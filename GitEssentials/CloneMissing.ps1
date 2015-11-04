##############################################
#	Title: Cloning All Lone non-Existing Locations (CALnEL)
#	Author: Kevin O'Neill
#	Date: July 5, 2012
##############################################

# Desc: Searches for present directories and clones the ones missing from RepoList.ini

# Use:	Execute from the GitEssentials folder.

##############################################

# Stores the git prompt as a vairable
$git = "C:\Git\cmd\git.cmd";

# Stores the list of repositories from RepoList.ini
$remoteRepos = Get-Content .\RepoList.ini;
.$git config core.autocrlf false;

# Goes up one directory level (GitEssentials folder -> CodeRepo folder)
cd ..;

# Searches for all folders. (eg: Project1, Project2, GitEssentials, ProjectZ)
$localRepos = ls -name -Exclude *.*;

# Iterates through all local repos
foreach ($repo in $remoteRepos)
{
	Write-Host "";	# Adding whitespace 
	
	# Compares local repos to the list of remote repos. If there is a match, do nothing, else clone the missin repo.
	if($localRepos -contains $repo)
	{
		Write-Host $repo "exists";
		if($repo -eq "GitEssentials")
		{
			cd .\$repo
			.$git branch --set-upstream master origin/master
			cd ..
		}
	}
	else
	{
		Write-Host "Cloning: " $repo;
		$remoteLocation = "http://CodeRepo@155.101.184.5/" + $repo + ".git"; 
		.$git clone $remoteLocation ; 
		Write-Host $repo "cloned.";
		
		if(Test-Path .\$repo -pathType Container)
		{
			cd .\$repo;
			.$git config core.autocrlf false;
			.$git checkout -b work_branch
			cd ..;
		}
	}
}

Write-Host "";	# Adding whitespace 

# Writes to command line that the script is done.
Write-Host "Done Cloning All Missing Repos.";

Write-Host "";	# Adding whitespace 

Write-host ""
Write-Host "Press any key to continue..."

# Waits for a key press to end program
$press = $host.UI.RawUI.ReadKey("NoEcho,IncludeKeyDown");

# EOF