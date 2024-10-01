## VIASH START
par_input=s3://openproblems-nextflow/work/f6/8565066aee4771cc2790b92b4ac660
par_aws_profile=op
par_output=debug
## VIASH END

aws s3 sync \
  ${par_aws_profile:+--profile $par_aws_profile} \
  "$par_input" "$par_output"

cd "$par_output"

sed -i 's#/home/ec2-user/miniconda/bin/aws#aws#g' .command.run

AWS_DEFAULT_PROFILE=${par_aws_profile} bash .command.run nxf_stage
