wget https://www.dropbox.com/sh/1xzt1rysy9v3a9l/AAAMlJBQG1OQFW4cjhp11Ex6a/LOW.tar.gz
wget https://bitbucket.org/grackle/grackle/src/default/input/CloudyData_noUVB.h5
tar xzvf LOW.tar.gz
mv LOW/*.dat ./
rmdir LOW
rm LOW.tar.gz