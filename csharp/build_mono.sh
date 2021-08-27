# !/usr/bin/bash
mcs -optimize+ main.cs Data.cs GILS_RVND.cs  && ./main.exe
rm *.exe
