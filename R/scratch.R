

################################################################################
#################   Fit VBC                     ################################
################################################################################

oc = cccma$rcm.c
mc = cccma$gcm.c
mp = cccma$gcm.p

m_control = list(mult = xmin = c(0, NaN, 0))
c_control = list(family_set = "onepar")
begin = Sys.time()
vbc = vine_correct(oc, mc, mp, margins_controls = m_control,
                   copula_controls = c_control) 
Sys.time() - begin
# debugonce(vine_correct)
#vbc = vine_correct(oc, mc, mp)
#vbc = vine_correct(oc, mc, mp, trunc0 = c(0.05, NA, 0))



escore(pseudo_obs(vbc), pseudo_obs(oc))
escore(pseudo_obs(mbcn.p), pseudo_obs(oc))
################################################################################
#################   Fit R2D2                    ################################
################################################################################
ubc = mapply(QDM, o.c = data.frame(oc), m.c = data.frame(mc), m.p = data.frame(mp),
       ratio=cccma$ratio.seq, trace=cccma$trace)
ubc = do.call(cbind, ubc["mhat.p",])
r2d2_BC = r2d2(oc, ubc)
r2d2_bc= r2d2_BC$r2d2_bc
colnames(r2d2_bc) = colnames(vbc)

################################################################################
#################   Plotting                    ################################
################################################################################

mbc_pair = grid.grabExpr(print(ggpairs(as.data.frame(mbcn.p),
        lower = list(continuous = wrap("density", alpha = 0.5, bins = 30)),
        upper = list(continuous = wrap("cor", method = "kendall")),
        title = "MBCn Correction")))
vbc_pair = grid.grabExpr(print(ggpairs(as.data.frame(vbc),
        lower = list(continuous = wrap("density", alpha = 0.5, bins = 30)),
        upper = list(continuous = wrap("cor", method = "kendall")),
        title = "Vine Correction")))
oc_pair = grid.grabExpr(print(ggpairs(as.data.frame(oc),
        lower = list(continuous = wrap("density", alpha = 0.5, bins = 30)),
        upper = list(continuous = wrap("cor", method = "kendall")),
        title = "Observed Data")))
r2_pair = grid.grabExpr(print(ggpairs(as.data.frame(r2d2_bc),
        lower = list(continuous = wrap("density", alpha = 0.5, bins = 30)),
        upper = list(continuous = wrap("cor", method = "kendall")),
        title = "R2D2 Correction")))

grid.arrange(oc_pair, vbc_pair)
grid.arrange(oc_pair, mbc_pair)
grid.arrange(oc_pair, r2_pair)
# ggpairs(as.data.frame(mp),
#         lower = list(continuous = wrap("density", alpha = 0.5, bins = 30)),
#         upper = list(continuous = wrap("cor", method = "kendall")),
#         title = "Model Data")


var = "tas"
par(mfrow=c(2,3))


acf(mbcn.p[, var])
acf(vbc[, var])
acf(r2d2_bc[, var])
acf(oc[, var])
# acf(mp[, var])
