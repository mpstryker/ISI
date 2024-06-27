# pgplot environment (csh)
setenv PGPLOT_DIR "/usr/local/pgplot/"
if ( $?LD_LIBRARY_PATH ) then 
    setenv LD_LIBRARY_PATH "${LD_LIBRARY_PATH}:/usr/local/lib/"
else
    setenv LD_LIBRARY_PATH /usr/local/lib/
endif
