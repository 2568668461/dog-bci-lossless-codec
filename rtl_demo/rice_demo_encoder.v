`timescale 1ns/1ps

// Minimal teaching core for the Zynq-7020 compression path.
// One frame is captured, then emitted as a bit stream:
//   first sample: 16-bit unsigned two's-complement, MSB first
//   remaining samples: zigzag(delta), Rice(K), MSB first
// The fixed K keeps this first demo deterministic and easy to verify.
module rice_demo_encoder #(
    parameter integer SAMPLE_COUNT = 16,
    parameter integer K = 2
) (
    input  wire        clk,
    input  wire        rst,
    input  wire        start,
    input  wire        in_valid,
    input  wire [15:0] in_data,
    output reg         in_ready,
    output reg         out_valid,
    output reg         out_bit,
    output reg         out_last,
    output reg         busy,
    output reg         done
);

    localparam [2:0] S_IDLE = 3'd0, S_CAPTURE = 3'd1, S_FIRST = 3'd2,
                     S_ZERO = 3'd3, S_ONE = 3'd4, S_REM = 3'd5,
                     S_DONE = 3'd6;

    reg [2:0] state;
    reg [15:0] samples [0:SAMPLE_COUNT-1];
    integer capture_idx;
    integer sample_idx;
    integer bit_idx;
    integer rem_idx;
    reg [16:0] quotient;
    reg [16:0] remainder;
    reg [15:0] previous;

    reg signed [17:0] delta_calc;
    reg [17:0] zigzag_calc;

    always @* begin
        in_ready = (state == S_CAPTURE);
        busy = (state != S_IDLE) && (state != S_DONE);
    end

    always @(posedge clk) begin
        if (rst) begin
            state <= S_IDLE;
            capture_idx <= 0;
            sample_idx <= 0;
            bit_idx <= 15;
            rem_idx <= 0;
            quotient <= 0;
            remainder <= 0;
            previous <= 0;
            out_valid <= 1'b0;
            out_bit <= 1'b0;
            out_last <= 1'b0;
            done <= 1'b0;
        end else begin
            out_valid <= 1'b0;
            out_last <= 1'b0;
            done <= 1'b0;

            case (state)
                S_IDLE: begin
                    if (start) begin
                        capture_idx <= 0;
                        state <= S_CAPTURE;
                    end
                end

                S_CAPTURE: begin
                    if (in_valid) begin
                        samples[capture_idx] <= in_data;
                        if (capture_idx == SAMPLE_COUNT-1) begin
                            bit_idx <= 15;
                            state <= S_FIRST;
                        end else begin
                            capture_idx <= capture_idx + 1;
                        end
                    end
                end

                S_FIRST: begin
                    out_valid <= 1'b1;
                    out_bit <= samples[0][bit_idx];
                    if (bit_idx == 0) begin
                        sample_idx <= 1;
                        previous <= samples[0];
                        if (SAMPLE_COUNT == 1) begin
                            state <= S_DONE;
                        end else begin
                            delta_calc = $signed({{2{samples[1][15]}}, samples[1]}) - $signed({{2{samples[0][15]}}, samples[0]});
                            zigzag_calc = ((delta_calc <<< 1) ^ (delta_calc >>> 17)) & 18'h1ffff;
                            quotient <= zigzag_calc >> K;
                            remainder <= zigzag_calc;
                            rem_idx <= K-1;
                            state <= S_ZERO;
                        end
                    end else begin
                        bit_idx <= bit_idx - 1;
                    end
                end

                S_ZERO: begin
                    if (quotient != 0) begin
                        out_valid <= 1'b1;
                        out_bit <= 1'b0;
                        quotient <= quotient - 1'b1;
                    end else begin
                        out_valid <= 1'b0;
                        state <= S_ONE;
                    end
                end

                S_ONE: begin
                    out_valid <= 1'b1;
                    out_bit <= 1'b1;
                    state <= (K == 0) ? S_DONE : S_REM;
                end

                S_REM: begin
                    out_valid <= 1'b1;
                    out_bit <= remainder[rem_idx];
                    if (rem_idx == 0) begin
                        previous <= samples[sample_idx];
                        if (sample_idx == SAMPLE_COUNT-1) begin
                            state <= S_DONE;
                        end else begin
                            delta_calc = $signed({{2{samples[sample_idx+1][15]}}, samples[sample_idx+1]}) - $signed({{2{samples[sample_idx][15]}}, samples[sample_idx]});
                            zigzag_calc = ((delta_calc <<< 1) ^ (delta_calc >>> 17)) & 18'h1ffff;
                            quotient <= zigzag_calc >> K;
                            remainder <= zigzag_calc;
                            rem_idx <= K-1;
                            sample_idx <= sample_idx + 1;
                            state <= S_ZERO;
                        end
                    end else begin
                        rem_idx <= rem_idx - 1;
                    end
                end

                S_DONE: begin
                    out_last <= 1'b1;
                    done <= 1'b1;
                    state <= S_IDLE;
                end

                default: state <= S_IDLE;
            endcase
        end
    end
endmodule
