####Running Polyploid Data in fastSTRUCTURE####

#all code in this script is performed on the HPC
#dataset has already been thinned

#convert data to fastSTRUCTURE compatible format
python3 Cochlearia_create_structure_file.py -v input_name/ -o output_name -s true
#input file should be in the format vcf - variant call format
#output file will be in the format str - structure format

###Run fastStructure with K values of 1 through 10 - do this 5 times each###
##run choose K each time. this is done to better understand stability##

#activate fastSTRUCTURE environment
conda activate /shared/conda/faststructure

#k1
python /path/to/faststructure/bin/structure.py -K 1 --input=/Path/To/File --output=Output_Name --format=str --full 

#k2
python /path/to/faststructure/bin/structure.py -K 2 --input=/Path/To/File --output=Output_Name --format=str --full 

#k3
python /path/to/faststructure/bin/structure.py -K 3 --input=/Path/To/File --output=Output_Name --format=str --full 

#k4
python /path/to/faststructure/bin/structure.py -K 4 --input=/Path/To/File --output=Output_Name --format=str --full 

#k5
python /path/to/faststructure/bin/structure.py -K 5 --input=/Path/To/File --output=Output_Name --format=str --full 

#k6
python /path/to/faststructure/bin/structure.py -K 6 --input=/Path/To/File --output=Output_Name --format=str --full 

#k7
python /path/to/faststructure/bin/structure.py -K 7 --input=/Path/To/File --output=Output_Name --format=str --full 

#k8
python /path/to/faststructure/bin/structure.py -K 8 --input=/Path/To/File --output=Output_Name --format=str --full 

#k9
python /path/to/faststructure/bin/structure.py -K 9 --input=/Path/To/File --output=Output_Name --format=str --full 

#k10
python /path/to/faststructure/bin/structure.py -K 10 --input=/Path/To/File --output=Output_Name --format=str --full 

#view - popfiles will be generated from the structures. these are text files

#k1
python /path/to/bin/distruct.py -K 1 --input=/Path/To/File --output=Path/To/Output/Location --title="K=1"  --popfile=/path/to/popfile_all.txt

#k2
python /path/to/bin/distruct.py -K 2 --input=/Path/To/File --output=Path/To/Output/Location --title="K=2"  --popfile=/path/to/popfile_all.txt

#k3
python /path/to/bin/distruct.py -K 3 --input=/Path/To/File --output=Path/To/Output/Location --title="K=3"  --popfile=/path/to/popfile_all.txt

#k4
python /path/to/bin/distruct.py -K 4 --input=/Path/To/File --output=Path/To/Output/Location --title="K=4"  --popfile=/path/to/popfile_all.txt

#k5
python /path/to/bin/distruct.py -K 5 --input=/Path/To/File --output=Path/To/Output/Location --title="K=5"  --popfile=/path/to/popfile_all.txt

#k6
python /path/to/bin/distruct.py -K 6 --input=/Path/To/File --output=Path/To/Output/Location --title="K=6"  --popfile=/path/to/popfile_all.txt

#k7
python /path/to/bin/distruct.py -K 7 --input=/Path/To/File --output=Path/To/Output/Location --title="K=7"  --popfile=/path/to/popfile_all.txt

#k8
python /path/to/bin/distruct.py -K 8 --input=/Path/To/File --output=Path/To/Output/Location --title="K=8"  --popfile=/path/to/popfile_all.txt

#k9
python /path/to/bin/distruct.py -K 9 --input=/Path/To/File --output=Path/To/Output/Location --title="K=9"  --popfile=/path/to/popfile_all.txt

#k10
python /path/to/bin/distruct.py -K 10 --input=/Path/To/File --output=Path/To/Output/Location --title="K=10"  --popfile=/path/to/popfile_all.txt


