##############################################
#	Title: Wacom Config Reset - Click (WaCoRe-Click)
#	Author: Kevin O'Neill
#	Date: July 23, 2012
##############################################

# Desc: Use to reset the Wacom configuration file.

# Use:	Run via powershell.

##############################################

# Gets the current script name
$scriptNameFull = ($MyInvocation.MyCommand).Name

# Empty string.
$emptyString = "";

# Truncates the full script name to exclude the extension
$scriptName = [regex]::replace($scriptNameFull, ".ps1", $emptyString);

$scriptName = $args[0];

# Initializes the various file and directory path needed to swap the two files.
$filePath   = ".\Wacom" + $scriptName + "\Wacom_Tablet.dat";
$targetDir  = "C:\Users\Administrator\AppData\Roaming\WTablet";
$targetFile = "C:\Users\Administrator\AppData\Roaming\WTablet\Wacom_Tablet.dat";

# Stops the Wacom Tablet Service
Stop-Service TabletServiceWacom

# Deletes the current configuration
Remove-Item $targetFile

# Copies in a known, desired, configuration
Copy-Item $filePath $targetDir

# Restarts the Wacom Tablet Service
Start-Service TabletServiceWacom

# EOF