library(data.table)

fb_traits <- readRDS("fb_data/fb_traits.rds")

# Split into trials 1 and 2
fb1 <- fb_traits[sdp != "progeny" | trial == 1]
fb2 <- fb_traits[sdp != "progeny" | trial == 2]

# Get list of parent IDs
parents <- fb_traits[sdp != "progeny", id]

# Split parents by appearance in each trial
parents1 <- parents[parents %in% fb1[, sort(unique(c(sire, dam)))]]
parents2 <- parents[parents %in% fb2[, sort(unique(c(sire, dam)))]]

# Clip out parents not present
fb1 <- fb1[id %in% parents1 | sdp == "progeny", ]
fb2 <- fb2[id %in% parents2 | sdp == "progeny", ]

# Remap IDss
fb1[, `:=`(sire = match(sire, unique(id)),
           dam = match(dam, unique(id)),
           id = match(id, unique(id)))]

fb2[, `:=`(sire = match(sire, unique(id)),
           dam = match(dam, unique(id)),
           id = match(id, unique(id)))]

# Fix groups for Trial 2
fb2[, group := group - 36L]

# Save
saveRDS(fb1, "fb_data/fb_traits1.rds")
saveRDS(fb2, "fb_data/fb_traits2.rds")
