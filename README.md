# Ilse's neuro repo

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

- change git's editor to be not vim:

```
git config --global core.editor "xdg-open"
```
