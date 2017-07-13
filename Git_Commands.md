# To Clone and keep updated 

1. Clone your fork:

git clone git@github.com:YOUR-USERNAME/YOUR-FORKED-REPO.git

2. Add remote from original repository in your forked repository:

cd into/cloned/fork-repo
git remote add upstream git://github.com/ORIGINAL-DEV-USERNAME/REPO-YOU-FORKED-FROM.git
git fetch upstream

3. Updating your fork from original repo to keep up with their changes:

git pull upstream master

# To pull and push to/from

git pull origin master



# To remove large files, GitHub suggests:

$ git rm --cached giant_file
#Stage our giant file for removal, but leave it on disk

git commit --amend -CHEAD
#Amend the previous commit with your change
#Simply making a new commit won't work, as you need
#to remove the file from the unpushed history as well
1
git push
#Push our rewritten, smaller commit


# Remove all local changes from your working copy, simply stash them:

git stash save --keep-index

If you don't need them anymore, you now can drop that stash:

git stash drop


