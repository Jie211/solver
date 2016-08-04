#!bin/zsh
for (( i=1; i<=4; i++))
do
  for (( j=10; j<=100; j+=10))
  do
    echo "===================start========================"
    cat ./CR/5000#30#1e-$i#$j#1/data.txt | grep bad
    cat ./CR/5000#30#1e-$i#$j#1/data.txt | grep ms
    echo "===================$i#$j over==================="
  done
done


