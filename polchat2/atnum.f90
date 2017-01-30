integer function atnum(atnam)

  use constants 

  character(len=typmax) :: atnam

 9000 format(' ERROR',/,&
             ' Cannot understand atom name ',(A))

  select case (trim(atnam))
    case ('H')  
      atnum =   1
    case ('He') 
      atnum =   2
    case ('Li') 
      atnum =   3
    case ('Be') 
      atnum =   4
    case ('B')  
      atnum =   5
    case ('C')  
      atnum =   6
    case ('N')  
      atnum =   7
    case ('O')  
      atnum =   8
    case ('F')  
      atnum =   9
    case ('Ne') 
      atnum =  10
    case ('Na') 
      atnum =  11
    case ('Mg') 
      atnum =  12
    case ('Al') 
      atnum =  13
    case ('Si') 
      atnum =  14
    case ('P')  
      atnum =  15
    case ('S')  
      atnum =  16
    case ('Cl') 
      atnum =  17
    case ('Ar') 
      atnum =  18
    case ('K')  
      atnum =  19
    case ('Ca') 
      atnum =  20
    case ('Sc') 
      atnum =  21
    case ('Ti') 
      atnum =  22
    case ('V')  
      atnum =  23
    case ('Cr') 
      atnum =  24
    case ('Mn') 
      atnum =  25
    case ('Fe') 
      atnum =  26
    case ('Co') 
      atnum =  27
    case ('Ni') 
      atnum =  28
    case ('Cu') 
      atnum =  29
    case ('Zn') 
      atnum =  30
    case ('Ga') 
      atnum =  31
    case ('Ge') 
      atnum =  32
    case ('As') 
      atnum =  33
    case ('Se') 
      atnum =  34
    case ('Br') 
      atnum =  35
    case ('Kr') 
      atnum =  36
    case ('Rb') 
      atnum =  37
    case ('Sr') 
      atnum =  38
    case ('Y')  
      atnum =  39
    case ('Zr') 
      atnum =  40
    case ('Nb') 
      atnum =  41
    case ('Mo') 
      atnum =  42
    case ('Tc') 
      atnum =  43
    case ('Ru') 
      atnum =  44
    case ('Rh') 
      atnum =  45
    case ('Pd') 
      atnum =  46
    case ('Ag') 
      atnum =  47
    case ('Cd') 
      atnum =  48
    case ('In') 
      atnum =  49
    case ('Sn') 
      atnum =  50
    case ('Sb') 
      atnum =  51
    case ('Te') 
      atnum =  52
    case ('I')  
      atnum =  53
    case ('Xe') 
      atnum =  54
    case ('Cs') 
      atnum =  55
    case ('Ba') 
      atnum =  56
    case ('La') 
      atnum =  57
    case ('Ce') 
      atnum =  58
    case ('Pr') 
      atnum =  59
    case ('Nd') 
      atnum =  60
    case ('Pm') 
      atnum =  61
    case ('Sm') 
      atnum =  62
    case ('Eu') 
      atnum =  63
    case ('Gd') 
      atnum =  64
    case ('Tb') 
      atnum =  65
    case ('Dy') 
      atnum =  66
    case ('Ho') 
      atnum =  67
    case ('Er') 
      atnum =  68
    case ('Tm') 
      atnum =  69
    case ('Yb') 
      atnum =  70
    case ('Lu') 
      atnum =  71
    case ('Hf') 
      atnum =  72
    case ('Ta') 
      atnum =  73
    case ('W')  
      atnum =  74
    case ('Re') 
      atnum =  75
    case ('Os') 
      atnum =  76
    case ('Ir') 
      atnum =  77
    case ('Pt') 
      atnum =  78
    case ('Au') 
      atnum =  79
    case ('Hg') 
      atnum =  80
    case ('Tl') 
      atnum =  81
    case ('Pb') 
      atnum =  82
    case ('Bi') 
      atnum =  83
    case ('Po') 
      atnum =  84
    case ('At') 
      atnum =  85
    case ('Rn') 
      atnum =  86
    case ('Fr') 
      atnum =  87
    case ('Ra') 
      atnum =  88
    case ('Ac') 
      atnum =  89
    case ('Th') 
      atnum =  90
    case ('Pa') 
      atnum =  91
    case ('U')  
      atnum =  92
    case ('Np') 
      atnum =  93
    case ('Pu') 
      atnum =  94
    case ('Am') 
      atnum =  95
    case ('Cm') 
      atnum =  96
    case ('Bk') 
      atnum =  97
    case ('Cf') 
      atnum =  98
    case ('Es') 
      atnum =  99
    case ('Fm') 
      atnum = 100
    case ('Md') 
      atnum = 101
    case ('No') 
      atnum = 102
    case ('Lr') 
      atnum = 103
    case ('Rf') 
      atnum = 104
    case ('Db') 
      atnum = 105
    case ('Sg') 
      atnum = 106
    case ('Bh') 
      atnum = 107
    case ('Hs') 
      atnum = 108
    case ('Mt') 
      atnum = 109
    case ('Ds') 
      atnum = 110
    case ('Rg') 
      atnum = 111
    case ('Cn') 
      atnum = 112
    case ('Nh') 
      atnum = 113
    case ('Fl') 
      atnum = 114
    case ('Mc') 
      atnum = 115
    case ('Lv') 
      atnum = 116
    case ('Ts') 
      atnum = 117
    case ('Og') 
      atnum = 118
    case default
      write(iout,9000) atnam
      stop
  end select

  return

end function
