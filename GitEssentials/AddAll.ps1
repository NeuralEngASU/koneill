##############################################
#	Title: Add All Changes (AAC)
#	Author: Kevin O'Neill
#	Date: July 26, 2012
##############################################

# Desc: Use to add all changes (additions, alterations, deletions) to the staged commit.

# Use:	Run inside the repo of choice

##############################################

# Stores the git prompt as a vairable
$git = "C:\Git\cmd\git.cmd";

# Adds all new and modified files to git gui
.$git add .

# Adds all tracked and deleted files to git gui
.$git add -u
