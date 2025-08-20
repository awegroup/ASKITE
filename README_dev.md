# Workflow
Following the git-branch-model: https://nvie.com/posts/a-successful-git-branching-model/
When merging develop to master, always ensure it is a stable version with a tag

## Implementing a new feature
### On GitHub
1. Create an issue on GitHub
2. OpenUp a branch from this issue and ensure it's made from the `develop` branch
### Locally
3. Update using the provided command (will include a `git fetch` and a `git checkout`)
4. Implement your new feature
### On GitHub
6. Use a pull request to merge, and inside of it add the current issue to close it -once merged
7. Once merged, delete this feature branch
8. Close issue
### Locally
9. Switch branch to develop, update branch information, pull in new info
```bash
git checkout develop
git fetch --prune
git pull origin develop
```
10. Delete branch
```bash
git branch -d <enter branch name>
```