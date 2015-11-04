##############################################
#	Title: Script for Updating Git Repos (SUGR)
#	Author: Kevin O'Neill
#	Date: May 19, 2012
##############################################

# Desc: Looks inside your Git repo and attempts to update existing repos.

# Use:	Execute from the GitEssentials folder.

##############################################

# Stores the git prompt as a vairable
$git = "C:\Git\cmd\git.cmd";

# Goes up one directory level (GitEssentials folder -> CodeRepo folder)
cd ..;

# Searches for all folders. (eg: Project1, Project2, GitEssentials, ProjectZ)
# Then for each folder, enters that folder and preforms a git pull
$localRepos = ls -name -Exclude *.* 

foreach ($repo in $localRepos)
{ 
	Write-Host "";	# Adding whitespace 
	
	Write-Host "Updating: " $repo;
	cd .\$repo; 	# Moves into folder
	.\PushRepo.ps1
	Write-Host "";
	Write-Host $repo "updated.";
	cd ..;	# Move up a directory
}

Write-Host "";	# Adding whitespace 

# Writes to command line that the script is done.
Write-Host "Done Updating All Existing Repos";

Write-Host "";	# Adding whitespace 

Write-Host "Press any key to continue..."
$press = $host.UI.RawUI.ReadKey("NoEcho,IncludeKeyDown");

# EOF
