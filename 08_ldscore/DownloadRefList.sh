# SFTP details
# This is the location your results will be uploaded
# Note: No individual data will be uploaded to this location during the course of the pipeline
# sftp_address="sscmv-filetran.epi.bris.ac.uk"
sftp_address="filetrn-scmv-d0.epi.bris.ac.uk"

if [ "$sftp_username" = "gh13047" ] || [ "$sftp_username" = "epzjlm" ] || [ "$sftp_username" = "godmc" ] || [ "$sftp_username" = "epwkb" ]
then
	#sftp_path="/srv/sftponly/GoDMC"
	sftp_path="/GoDMC"
else
	sftp_path="/GoDMC"
fi


sftp_username="ehannon"


sftp ${sftp_username}@${sftp_address}:${sftp_path}/resources/phase2 <<EOF
get lists_17.tgz
get lists_17.tgz.md5sum
EOF

echo "Checking download integrity"

md5sum -c lists_17.tgz.md5sum

echo "Extracting"

mkdir resources
tar xzf lists_17.tgz -C resources
#rm lists_17.tgz*

