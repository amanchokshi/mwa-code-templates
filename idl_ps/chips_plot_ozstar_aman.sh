echo 'starting'
echo .compile \'/fred/oz048/achokshi/mwa_dipole/simulations/full_sim_II/ps/plot_chipsout_aman.pro\' > idl_batch_aman
echo FILE_MKDIR, \'${1}\' >> idl_batch_aman
echo CD, \'${1}\' >> idl_batch_aman 
echo .go >> idl_batch_aman
echo xx_0.iter.${1} >> idl_batch_aman
echo yy_0.iter.${1} >> idl_batch_aman
echo 1 >> idl_batch_aman
echo $1 >> idl_batch_aman

# Creates new script to rotate idl png images by 180 deg
cat > rotate$1.sh <<EOF
for file in ./$1/*.png; do
    convert "\$file" -rotate 180 "\${file%.png}".png
done
EOF
chmod +x rotate$1.sh

echo 'done'
