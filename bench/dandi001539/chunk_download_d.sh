#!/bin/bash
# 分段下载大鼠皮层 NWB（在 D 盘运行，绕开单连接限时 + C 盘满）
set -e
cd /d/bcidata
TARGET=3374757242
FILE=Rat06_Insertion2_Depth1.nwb
CHUNK=157286400  # 150 MB (~40s @4MB/s)

curl -s -o tmp_head.bin -w "%{redirect_url}" "https://huggingface.co/datasets/rokaijano/rat_cortical_128ch/resolve/main/raw_data/Rat06_Insertion2_Depth1.nwb" > cdn_url.txt
CDN=$(cat cdn_url.txt)

get_size() {
  ls -la "$FILE" 2>/dev/null | awk '{print $5}' || echo 0
}

POS=$(get_size)
echo "start from $POS / $TARGET"
while [ "$POS" -lt "$TARGET" ]; do
  END=$((POS + CHUNK - 1))
  if [ "$END" -ge "$TARGET" ]; then END=$((TARGET - 1)); fi
  for attempt in 1 2 3 4 5 6 7 8; do
    curl -s -r "$POS-$END" -o chunk_tmp.bin "$CDN" && break
    echo "retry $attempt for range $POS-$END"
    curl -s -o tmp_head.bin -w "%{redirect_url}" "https://huggingface.co/datasets/rokaijano/rat_cortical_128ch/resolve/main/raw_data/Rat06_Insertion2_Depth1.nwb" > cdn_url.txt
    CDN=$(cat cdn_url.txt)
    "C:/Users/MI/AppData/Local/Programs/Python/Python312/python.exe" -c "import time;time.sleep(3)"
  done
  GOT=$(ls -la chunk_tmp.bin | awk '{print $5}')
  EXPECT=$((END - POS + 1))
  if [ "$GOT" != "$EXPECT" ]; then
    echo "size mismatch: got $GOT expect $EXPECT, aborting at POS=$POS"
    exit 1
  fi
  cat chunk_tmp.bin >> "$FILE"
  POS=$(get_size)
  echo "progress: $POS / $TARGET ($(( POS * 100 / TARGET ))%)"
done
echo "DOWNLOAD_COMPLETE $POS bytes"
