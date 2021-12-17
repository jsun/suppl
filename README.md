# suppl


```bash
git clone suppl
git clone tiramisu

cd suppl
mkdir dragonfly
touch dragonfly/.gitkeep

git add -A tiramisu
git commit -m 'prepare to merge tiramisu.'
git remote add tiramisu ~/tiraimsu
git fetch tiramisu
git merge --allow-unrelated-histories -X subtree=tiraimsu tiramisu/master

git mv dragonfly 10.3389/fevo.2021.762173
rm 10.3389/fevo.2021.762173/.gitkeep
git rm 10.3389/fevo.2021.762173/.gitkeep
git commit -m 'change dirname to doi'

git log --graph --all
```


```bash
IFS=$'\n';

objects=`git verify-pack -v .git/objects/pack/pack-*.idx | grep -v chain | sort -k3nr | head -n 100`

output="size,pack,SHA,location"
for y in $objects
do
    # extract the size in bytes
    size=$((`echo $y | cut -f 5 -d ' '`/1024))
    # extract the compressed size in bytes
    compressedSize=$((`echo $y | cut -f 6 -d ' '`/1024))
    # extract the SHA
    sha=`echo $y | cut -f 1 -d ' '`
    # find the objects location in the repository tree
    other=`git rev-list --all --objects | grep $sha`
    #lineBreak=`echo -e "\n"`
    output="${output}\n${size},${compressedSize},${other}"
done

echo -e $output | column -t -s ', '
```


```bash
git filter-branch -f --tree-filter "rm -rf tmp" -- --all
git checkout master && git push -f origin master
```

