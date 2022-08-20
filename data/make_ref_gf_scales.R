ref_gf_scales <- rbind(
  data.frame(val=0.1465, dval=0.0024,
             label="$(t_0^{\\mathrm{plaq}})^{1/2}/a$", 
             obs="tsqEplaq",
             obslabel="$\\langle t^2 E_{\\mathrm{plaq}}(t) \\rangle$",
             Nf="$Nf=2+1$ (Wilson)",
             ref="Borsanyi et al., JHEP 09 (2012) 010, High-precision scale setting in lattice QCD"),
  data.frame(val=0.1465, dval=0.0024,
             label="$(t_0^{\\mathrm{plaq}})^{1/2}/a$", 
             obs="tsqEsym",
             obslabel="$\\langle t^2 E_{\\mathrm{sym}}(t) \\rangle$",
             Nf="$Nf=2+1$ (Wilson)",
             ref="Borsanyi et al., JHEP 09 (2012) 010, High-precision scale setting in lattice QCD"),
  data.frame(val=0.1755, dval=0.019,
             label="(w_0^{\\mathrm{sym}}/a)^2",
             obs="Wsym",
             obslabel="$\\langle t \\frac{d}{dt} [ t^2 E(t) ] \\rangle$",
             Nf="$Nf=2+1$ (Wilson)",
             ref="Borsanyi et al., JHEP 09 (2012) 010, High-precision scale setting in lattice QCD")
)
save(ref_gf_scales, file = "ref_gf_scales.RData")
