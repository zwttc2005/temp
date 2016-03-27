    //-------------------------------------------------------------------
    // Load the refprop thermodynamic libaray 
    //-------------------------------------------------------------------
    std::string err;
    bool loaded_REFPROP = load_REFPROP(err);
    long ierr = 0, nc = 1; 
    char herr[255], hfld[] = "CO2.FLD", hhmx[] = "HMX.BNC", href[] = "DEF";

    SETUPdll(nc,hfld,hhmx,href,ierr,herr,10000,255,3,255);

    Info<< "error massage:" << herr << nl << endl; 
    //-------------------------------------------------------------------	
