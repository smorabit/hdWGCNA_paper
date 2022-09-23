#!/bin/bash

# inline arguments:
indir=$1
ref=$2
cores=$3
mem=$4

# make a tmp directory:
cd $indir
mkdir tmp/

# what is the name of the whitelist file present in the current dir?
whitelist=$(ls ${indir}/*whitelist*)
echo $whitelist

################################################################################
# bustools on the output.bus file with spliced + unspliced
#
# output.bus ->
#
# bustools sort
# bustools inspect
# bustools correct
# bustools sort
#
# -> output.unfiltered.bus
################################################################################

echo "Running commands on output.bus"

# sort spliced
bustools sort \
  -t $cores \
  -m $mem \
  -o  tmp/output.s.bus \
  output.bus

# inspect spliced
bustools inspect \
  -o tmp/inspect.json \
  -e matrix.ec \
  -w $whitelist \
  tmp/output.s.bus

bustools correct \
  -o tmp/output.s.c.bus \
  -w $whitelist \
  tmp/output.s.bus

# sort again
bustools sort \
  -t $cores \
  -m $mem \
  -o  tmp/output.unfiltered.bus \
  tmp/output.s.c.bus



################################################################################
# get the spliced matrix
#
# output.unfiltered.bus ->
#
# bustools capture
# bustools sort
# bustools inspect
# bustools count
#
# -> spliced.mtx
################################################################################

echo "Processing spliced data"

# bustools capture for spliced
bustools capture \
  -s \
  -e matrix.ec \
  -t transcripts.txt \
  --output tmp/spliced.bus \
  --capture $ref/intron_t2c.txt \
  tmp/output.unfiltered.bus # input file

# sort spliced
bustools sort \
  -t $cores \
  -m $mem \
  -o  tmp/spliced.unfiltered.bus \
  tmp/spliced.bus

# inspect spliced
bustools inspect \
  -o tmp/inspect.spliced.json \
  -e matrix.ec \
  -w $whitelist \
  tmp/spliced.unfiltered.bus

# count spliced
bustools count \
  -o tmp/spliced \
  -g $ref/t2g.txt \
  -e matrix.ec \
  -t transcripts.txt \
  --genecounts \
  tmp/spliced.unfiltered.bus


################################################################################
# get the unspliced matrix
#
# output.unfiltered.bus ->
#
# bustools capture
# bustools sort
# bustools inspect
# bustools count
#
# -> unspliced.mtx
################################################################################

echo "Processing unspliced data"


# bustools capture for unspliced
bustools capture \
  -s \
  -e matrix.ec \
  -t transcripts.txt \
  --output tmp/unspliced.bus \
  --capture $ref/cdna_t2c.txt \
  tmp/output.unfiltered.bus # input file

# sort unspliced
bustools sort \
  -t $cores \
  -m $mem \
  -o  tmp/unspliced.unfiltered.bus \
  tmp/unspliced.bus

# inspect unspliced
bustools inspect \
-o tmp/inspect.unspliced.json \
  -e matrix.ec \
  -w $whitelist \
  tmp/unspliced.unfiltered.bus

# count unspliced
bustools count \
  -o tmp/unspliced \
  -g $ref/t2g.txt \
  -e matrix.ec \
  -t transcripts.txt \
  --genecounts \
  tmp/unspliced.unfiltered.bus
