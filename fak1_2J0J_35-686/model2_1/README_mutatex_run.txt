#source virtual environment
. /usr/local/envs/mutatex/bin/activate
#replace file.pdb with the file name of your pdb
tsp -N 4 nohup mutatex file.pdb -m mutation_list.txt --foldx-version=suite5 --np 4 -x /usr/local/foldx5_2020/foldx -B -C deep  -v -L -c -l >& log &
