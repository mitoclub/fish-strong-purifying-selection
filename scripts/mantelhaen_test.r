# CMH test for all families using weight values (non-significant)

exp5_lab <-
  array(c(138, 91, 92, 22,        #1x1
          117, 105, 75, 51,       #1x2
          132, 121, 118, 90,      #2x1
          149, 111, 130, 135),    #2x2
        dim = c(2, 2, 4),
        dimnames = list(
          Mutagen = c("M-", "M+"),
          Heat_shock = c("H-", "H+"),
          Family = c("1", "2", "3", "4")))

mantelhaen.test(exp5_lab)
