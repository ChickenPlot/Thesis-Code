#!/bin/bash

# sudo apt install pigz
# activate conda 

# Set the directory containing your files
data_dir=../X204SC22090037-Z01-F002

rawdata_dir=../X204SC22090037-Z01-F002/01.RawData
out_dir=../output


mkdir ../output/assembly_genomeref
assembly_dir=${out_dir}/assembly_genomeref

refs_dir=../../arabidopsis_reference/transcriptome

# Use an array to store the filenames
datafiles=($(ls ${rawdata_dir} -1v))

# Loop through the files two by two
for ((i=0; i<${#datafiles[@]}; i+=2)); do

    fwd_file="${datafiles[i]}" # forward strand
    rv_file="${datafiles[i + 1]}" # reverse strand

    echo "Start unzipping files"
    
    echo "Forward: $fwd_file"
    pigz -d ${rawdata_dir}/${fwd_file}  # unzip
    fwd_file_unzipped=${fwd_file%.gz} # assign a name without the .gz
     
    echo "Reverse: $rv_file"
    pigz -d ${rawdata_dir}/${rv_file}
    rv_file_unzipped=${rv_file%.gz}

    out=${assembly_dir}/$(basename ${fwd_file%%_1.fq.gz})
    
    echo
    # STAR 
    echo "Running STAR"
    STAR --runMode alignReads --genomeLoad  NoSharedMemory --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --genomeDir "${refs_dir}/indexed_genome" --readFilesIn ${rawdata_dir}/${fwd_file_unzipped} ${rawdata_dir}/${rv_file_unzipped} --runThreadN 8 --outFileNamePrefix ${out}

    echo
    #rizippare
    echo "Zipping files back"
    echo "Forward"
    pigz ${rawdata_dir}/${fwd_file_unzipped}

    echo "Reverse"
    pigz ${rawdata_dir}/${rv_file_unzipped}
    
    # Add a blank line for separation
   echo
done

echo "Done processing assembly!"
