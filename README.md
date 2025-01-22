# Ilse's neuro repo

## Most used commands

- download/sync new repo: `git clone <git@XXX>`
- add file/file change to be tracked by git: `git add <filename>`
- commit currently added set of changes: `git commit -m "brief explanation"`
- upload new commits: `git push`
- download new commits: `git pull`
- check for new commits without downloading them: `git fetch`

## Git configuration

- Add user name and mail:

```
git --global config user.email 125893237+iklinkhamer@users.noreply.github.com
git --global config user.name ilse
```

- Create ssh key: `ssh-keygen`
	- Subsequently must be uploaded as authentication ssh key to github (or similar)
	
- Start and upload a new repository
	- start repo locally: `git init`
	- configure remote to upload to: `git remote add <git@github.com:XXXX>`
	- make first push: `git push -u origin master/Main`

- change git's editor to be not vim (requires nano to be installed):

```
git config --global core.editor "nano"
```
