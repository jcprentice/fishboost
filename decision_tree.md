if Tsym == "no" (=> Trec == "no")
    if parasites == FALSE
        if first donor in group 71 (id=1805)
            [E:0],[E|I:Tmax]
        if any other donor or recipient
            [S:0],[S:Tmax]
    if parasites == TRUE
        if donor
            [E:0],[E|I:Tmax]
        if recipient
            [S:0],[E|I:1],[E|I:Tmax]
if Tsym == "."
    if donor
       [E:0],[E|I|D:1],[E|I|D:Trec-1],[R:Trec]
    if recipient
       [S:0],[E|I|D:1],[E|I|D:Trec-1],[R:Trec]
if Tsym == a number
    if donor
       [E:0],[E|I:1],[E|I:Tsym-1],[D:Tsym],[D:Trec-1],[R:Trec],[R:Tmax]
    if recipient
       [S:0],[E|I:1],[E|I:Tsym-1],[D:Tsym],[D:Trec-1],[R:Trec],[R:Tmax]

