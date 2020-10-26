setwd('/cluster/projects/hansengroup/sujunc/scRNA/primary/scran/revision/');
#data from /Users/sujunchen/Dropbox/me/singleCell/figures/RF/elements/KLK3_celline_coculture_R_ttest.pzfx
# celline EV KLK3 CellineExo
a1 <- c(0.866708, 1.153792);
a2 <- c(562.5307, 468.6868, 634.7763);
t.test(a1, a2)$p.value
# cell line KLK3
a1 <- c(0.92649029, 1.07934213);
a2 <- c(3768.40755, 4762.18733);
t.test(a1, a2)$p.value;
# co-culture with Exo
a1 <- c(0.918785, 1.170684, 0.940761);
a2 <- c(146.6161, 119.5005, 128.3762);
t.test(a1, a2)$p.value;
#data from /Users/sujunchen/Dropbox/me/singleCell/figures/F3/elements/KLK3_PMSAsort.pzfx
# PMSA +/- KLK3
a1 <- rep(1, 3);
a2 <- c(71.42, 88.02, 84.73);
t.test(a1, a2)$p.value;
#cell line co-culture 
a1 <- rep(1, 3);
a2 <- c(29097.41, 49219.49, 52681);
t.test(a1, a2)$p.value;
