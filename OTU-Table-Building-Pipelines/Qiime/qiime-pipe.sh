module load qiime
prefix="../result_C101HW18071131/01.CleanData/"
for file in *;
do
    #echo $file
    foo=${file#"$prefix"}
    echo "$foo";
    bash silva-qiime-pipe.sh $file/$foo.fna $foo
done

