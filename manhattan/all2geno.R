allel2geno <- function(x){
  
  if (x == "A"){
    return("A/A")
    
  }else if (x == "G"){
    return("G/G")
    
  }else if (x == "C"){
    return("C/C")
    
  }else if (x == "T"){
    return("T/T")
    
  }else if (x == "R"){
    return("A/G")
    
  }else if (x == "Y"){
    return("C/T")
    
  }else if (x == "S"){
    return("G/C")
    
  }else if (x == "W"){
    return("A/T")
    
  }else if (x == "K"){
    return("G/T")
    
  }else if (x == "M"){
    return("A/C")
    
  }else if (x == "B"){
    return("C/G/T")
    
  }else if (x == "D"){
    return("A/G/T")
    
  }else if (x == "H"){
    return("A/C/T")
    
  }else if (x == "V"){
    return("A/C/G")
    
  }else if (x == "N"){
    return("N")
    
  } else if (x == "-" || x == ".") {
    return("gap")}
  
}
