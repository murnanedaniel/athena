rm log
for file in *.dat; do
    echo -e comparing "$file" and "$(basename ${file%.*}).txt"
    diff "$file" "$(basename ${file%.*}).txt" > log
done
echo -e results stored in log file
