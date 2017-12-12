
zcat ../results/16/16*.txt.gz |grep cg03636183 >cg03636183.txt
zcat ../results/16/16_1.txt.gz |head -n1 >header.txt
cat cg03636183.txt >>header.txt
mv header.txt cg03636183.txt
