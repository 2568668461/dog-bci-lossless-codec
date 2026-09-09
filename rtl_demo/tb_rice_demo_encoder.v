`timescale 1ns/1ps

module tb_rice_demo_encoder;
    localparam integer N = 8;
    reg clk = 1'b0;
    reg rst = 1'b1;
    reg start = 1'b0;
    reg in_valid = 1'b0;
    reg [15:0] in_data = 16'd0;
    wire in_ready, out_valid, out_bit, out_last, busy, done;
    integer i, out_file;
    reg signed [15:0] samples [0:N-1];

    rice_demo_encoder #(.SAMPLE_COUNT(N), .K(2)) dut (
        .clk(clk), .rst(rst), .start(start),
        .in_valid(in_valid), .in_data(in_data), .in_ready(in_ready),
        .out_valid(out_valid), .out_bit(out_bit), .out_last(out_last),
        .busy(busy), .done(done)
    );

    always #5 clk = ~clk;

    always @(posedge clk)
        if (out_valid) $fwrite(out_file, "%0d\n", out_bit);

    initial begin
        samples[0] = 16'sd100; samples[1] = 16'sd101;
        samples[2] = 16'sd99;  samples[3] = 16'sd99;
        samples[4] = 16'sd104; samples[5] = 16'sd100;
        samples[6] = 16'sd100; samples[7] = 16'sd101;
        out_file = $fopen("rtl_bits.txt", "w");
        repeat (2) @(posedge clk);
        rst <= 1'b0;
        @(posedge clk); start <= 1'b1;
        @(posedge clk); start <= 1'b0;
        for (i = 0; i < N; i = i + 1) begin
            @(posedge clk);
            in_valid <= 1'b1;
            in_data <= samples[i];
            while (!in_ready) @(posedge clk);
        end
        @(posedge clk); in_valid <= 1'b0; in_data <= 16'd0;
        wait (done);
        @(posedge clk);
        $fclose(out_file);
        $display("PASS: RTL stream written to rtl_bits.txt");
        $finish;
    end
endmodule
