nwjs_dir="d:/Programming/nwjs-v0.12.3-win-x64/nwjs-v0.12.3-win-x64"
build_app_dir="./app_build"
app_name="GlycReSoft"

echo $nwjs_dir
echo $build_app_dir

rm -rf $build_app_dir
mkdir -p $build_app_dir

# zip all files to nw archive
zip -r $app_name.nw ./package.json ./tk_local_interface.py ./static/nwjs ./static/dist

# copy nw.pak from current build node-webkit
cp $nwjs_dir/nw.pak ./nw.pak
cp $nwjs_dir/icudtl.dat ./icudtl.dat

# compilation to executable form
cat $nwjs_dir/nw.exe ./$app_name.nw > $build_app_dir/$app_name.exe && chmod +x $build_app_dir/$app_name.exe

# move nw.pak to build folder
mv ./nw.pak $build_app_dir/nw.pak
mv ./icudtl.dat $build_app_dir/icudtl.dat
# cp -r ./static/dist/* $build_app_dir

# remove $app_name.nw
rm ./$app_name.nw

# run application
$build_app_dir/$app_name.exe
