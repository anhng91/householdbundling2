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



while [ ! -f ../../householdbundling_estimate/estimate_1723941858.rds ]
do
    sleep 0.1m
done

for i in {1..10}
do
   sudo Rscript estimation/estimation.R
done

scp -i household-bundling-96.pem ec2-user@ec2-18-220-216-157.us-east-2.compute.amazonaws.com:/home/ec2-user/householdbundling_estimate/estimate_1723885745.rds .

sudo R -e "install.packages('tidyverse', repos='http://cran.rstudio.com/')"
sudo R -e "install.packages('pracma', repos='http://cran.rstudio.com/')"
sudo R -e "install.packages('randtoolbox', repos='http://cran.rstudio.com/')"
sudo R -e "install.packages('mvtnorm', repos='http://cran.rstudio.com/')"
sudo R -e "install.packages('ggplot2', repos='http://cran.rstudio.com/')"
sudo R -e "install.packages('lfe', repos='http://cran.rstudio.com/')"
sudo R -e "install.packages('plm', repos='http://cran.rstudio.com/')"
sudo R -e "install.packages('stargazer', repos='http://cran.rstudio.com/')"
sudo R -e "install.packages('randomForest', repos='http://cran.rstudio.com/')"
sudo R -e "install.packages('splitfngr', repos='http://cran.rstudio.com/')"
sudo R -e "install.packages('Hmisc', repos='http://cran.rstudio.com/')"
sudo R -e "install.packages('gridExtra', repos='http://cran.rstudio.com/')"
sudo R -e "install.packages('fastDummies', repos='http://cran.rstudio.com/')"


#!/bin/bash
#install R
sudo yum install -y R

#install RStudio-Server 1.0.153 (2017-07-20)
wget https://download2.rstudio.org/rstudio-server-rhel-1.0.153-x86_64.rpm
sudo yum install -y --nogpgcheck rstudio-server-rhel-1.0.153-x86_64.rpm
rm rstudio-server-rhel-1.0.153-x86_64.rpm

#install shiny and shiny-server (2017-08-25)
R -e "install.packages('shiny', repos='http://cran.rstudio.com/')"
wget https://download3.rstudio.org/centos5.9/x86_64/shiny-server-1.5.4.869-rh5-x86_64.rpm
sudo yum install -y --nogpgcheck shiny-server-1.5.4.869-rh5-x86_64.rpm
rm shiny-server-1.5.4.869-rh5-x86_64.rpm

sudo yum install -y freetype-devel libpng-devel libtiff-devel libjpeg-devel
sudo yum install -y libcurl-devel
sudo yum install -y openssl-devel
sudo yum install -y harfbuzz-devel fribidi-devel
sudo yum install -y tmux
sudo yum install -y git 

#install necessary packages 
sudo R -e "install.packages('tidyverse', repos='http://cran.rstudio.com/')"
sudo R -e "install.packages('pracma', repos='http://cran.rstudio.com/')"
sudo R -e "install.packages('randtoolbox', repos='http://cran.rstudio.com/')"
sudo R -e "install.packages('mvtnorm', repos='http://cran.rstudio.com/')"
sudo R -e "install.packages('ggplot2', repos='http://cran.rstudio.com/')"
sudo R -e "install.packages('lfe', repos='http://cran.rstudio.com/')"
sudo R -e "install.packages('plm', repos='http://cran.rstudio.com/')"
sudo R -e "install.packages('stargazer', repos='http://cran.rstudio.com/')"
sudo R -e "install.packages('randomForest', repos='http://cran.rstudio.com/')"
sudo R -e "install.packages('splitfngr', repos='http://cran.rstudio.com/')"
sudo R -e "install.packages('Hmisc', repos='http://cran.rstudio.com/')"
sudo R -e "install.packages('gridExtra', repos='http://cran.rstudio.com/')"
sudo R -e "install.packages('fastDummies', repos='http://cran.rstudio.com/')"
sudo R -e "install.packages('devtools', repos='http://cran.rstudio.com/')"


mkdir personal_data
cd personal_data

sudo git clone https://github.com/anhng91/householdbundling2.git 

cd householdbundling2

sudo R -e "devtools::install(upgrade='never')"












scp -i household-bundling.pem -r ec2-user@ec2-52-202-142-124.compute-1.amazonaws.com:/home/ec2-user/householdbundling_estimate/. ./householdbundling_estimate/.
scp -i household-bundling.pem -r ec2-user@ec2-52-202-142-124.compute-1.amazonaws.com:/home/ec2-user/Obj_for_manuscript/. ./Obj_for_manuscript/.


declare -a arr=("ec2-18-217-219-213.us-east-2.compute.amazonaws.com" "ec2-3-145-5-212.us-east-2.compute.amazonaws.com" "ec2-13-58-54-160.us-east-2.compute.amazonaws.com" "ec2-3-142-77-92.us-east-2.compute.amazonaws.com" "ec2-3-143-253-238.us-east-2.compute.amazonaws.com")
for ip in "${arr[@]}"
do
scp -i household-bundling-96.pem -r ec2-user@$ip:/home/ec2-user/householdbundling_estimate/. ./householdbundling_estimate/.
scp -i household-bundling-96.pem -r ec2-user@$ip:/home/ec2-user/Obj_for_manuscript/. ./Obj_for_manuscript/.            
done





