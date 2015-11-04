% Compare the Sync signals from XlTec and the modified TDT signal

timeMatXlTec = timeMat(2:end-1);
timeMatTDT = timeStampTDT(2:end-1);

idx = strfind(timeMatXlTec, timeMatTDT);