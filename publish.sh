if [ $# -eq 0 ]
then
	echo "No version bump. To bump version, pass major/minor/patch "
else
	bumpversion --allow-dirty $1
fi

rm -rf build
rm -rf dist
rm -rf *.egg-info

python3 setup.py sdist
python3 setup.py bdist_wheel --universal

read -rp "Upload? [yes/no] "

if [[ ${REPLY} =~ ^(yes)$ ]]; then
	twine upload dist/*
else
	echo "Upload canceled "
fi
