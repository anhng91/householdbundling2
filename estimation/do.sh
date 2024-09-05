declare -a list_servers=("ec2-18-216-94-114" "ec2-18-222-86-243" "ec2-3-129-13-57" "ec2-3-144-155-200" "ec2-3-23-127-53")


for server in "${list_servers[@]}"
do
scp -i ../../household-bundling-96.pem -r  ec2-user@"$server".us-east-2.compute.amazonaws.com:/home/ec2-user/householdbundling_estimate/. ../../householdbundling_estimate/
done 

scp -i ../../household-bundling.pem -r  ec2-user@ec2-54-197-6-122.compute-1.amazonaws.com:/home/ec2-user/householdbundling_estimate/. ../../householdbundling_estimate/



sudo rm -r householdbundling_estimate
cd personal_data 
sudo git clone https://github.com/anhng91/householdbundling2.git 
cd householdbundling2
tmux new -t mysession 
sudo Rscript estimation/estimation 