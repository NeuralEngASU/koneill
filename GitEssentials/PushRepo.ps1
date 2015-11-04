##############################################
#	Title: Pull/Rebase/Merge/Push (PuReMerPu)
#	Author: Kevin O'Neill
#	Date: July 23, 2012
##############################################

# Desc: Use to push and update current repo

# Use:	Run inside the repo of choice

##############################################

# Stores the git prompt as a vairable
$git = "C:\Git\cmd\git.cmd";

# Get the current path and convert to a string.
$currentPath = Get-Location;
$currentPath = [string]$currentPath;

# Go up one directory level
cd ..

# Get the parent directory path and convert to a string. Concantenates a '\' to the end of the string.
$parentPath = Get-Location;
$parentPath = [string]$parentPath + "\";

# Empty string.
$emptyString = "";

# Deletes the parent directory from the current directory in order to isolate the current folder name.
$myFolder = [regex]::replace($currentPath, [regex]::Escape($parentPath), $emptyString);

cd .\$myFolder

Write-Host "";
Write-Host "Updating: " $myFolder;

# Checkout to master and update master.
.$git checkout master
.$git pull

# Checkout to work_branch and rebase work_branch
.$git checkout work_branch
.$git rebase master

# Checkout to master and attempt to merge in work_branch
.$git checkout master
.$git merge work_branch

# Deletes the branch, pushes changes to the remote repo, and creates a new branch.
.$git branch -d work_branch
.$git push origin master
.$git checkout -b work_branch

Write-Host "";
Write-Host "Done pushing: " $myFolder;

Write-host ""
Write-Host "Press any key to continue..."

# Waits for a key press to end program
$press = $host.UI.RawUI.ReadKey("NoEcho,IncludeKeyDown");

# EOF