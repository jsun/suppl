library(rgl)



fugu.id <- '3'


## Argument Settings
argv <- commandArgs(TRUE)
dat.file    <- paste0('fugustats_cm.Fugu', fugu.id, '.vdj.freq.tsv')
sample.name <- paste0('Fugu ', fugu.id)
prob        <- TRUE
fig.name    <- paste0('Fig_Cm_VDJ3D.Fugu', fugu.id,'.version2')



## Functions
.colramp <- function(val, prob) {
    cols <- c("#deebf7", "#08306b")
    val <- val / 15
    val[val > 1] <- 1  # round 15.060695 to 15, because we will set 15 as a maximum
    x <- colorRamp(cols)(val)
    rgb(x[, 1], x[, 2], x[, 3], maxColorValue = 255)
}

vlevels <- c("v1.1", "v1.2", "v1.3", "v1.4", "v1.7", "v1.8",
             "v1.12", "v1.13", "v1.14", "v1.15", "v1.16", "v1.17", "v1.18", "v1.21",
             "v2.1", "v2.2", "v2.3", "v2.6", "v2.7", "v2.8", "v2.9",
             "v2.10", "v2.11", "v2.12", "v2.15", "v2.16", "v2.17", "v2.18", "v2.19", "v2.20",
             "v3.2", "v3.3")

if(F) {
vlevels <- c("v2.1", "v1.1", "v2.2", "v1.2", "v1.3",
             "v2.3", "v1.4", "v3.2", "v1.7", "v3.3",
             "v2.6", "v1.8", "v2.7", "v2.8", "v2.9", "v2.10",
             "v2.11", "v1.12", "v2.12", "v1.13", "v1.14", "v1.15",
             "v2.15", "v1.16", "v2.16", "v2.17", "v1.17", "v2.18",
             "v1.18", "v2.19", "v2.20", "v1.21")
}

jlevels <- c("Jm1", "Jm2", "Jm3", "Jm4", "Jm5")
dlevels <- c("Dm1", "Dm2", "Dm3", "Dm4", "Dm5", "Dm6", "Dm?")

label2axis <- function(a, g) {
    if (g == 'v') f <- vlevels
    if (g == 'd') f <- dlevels
    if (g == 'j') f <- jlevels
    a.factor <- factor(as.character(a), levels = f)
    as.numeric(a.factor)
}


## Read Data
dat <- read.table(dat.file, skip = 1)
dat <- dat[, -1]
dat[, 2] <- as.character(dat[, 2])
dat[is.na(dat[, 2]), 2] <- "Dm?"
dat <- dat[!is.na(dat[, 2]), ]
colnames(dat) <- c("V", "D", "J", "value")
if (prob == "TRUE") dat[, 4] <- 100 * dat[, 4] / sum(dat[, 4])


## 3D Plot
if (prob == "FALSE") {
    dat$value_adj <- dat$value ^ (1/4) / 10
} else {
    dat$value_adj <- (dat$value * 1) ^ (1/4) / 1
}
## add dummy data to adjust axes
dat <- rbind(dat, data.frame(V = 'v3.3', D = 'Dm6', J = 'Jm5', value = 0, value_adj = 0))
cols <- .colramp(dat$value, prob)

aspect3d(1, 1, 1)
plot3d(label2axis(dat$V, "v"), label2axis(dat$D, "d"), label2axis(dat$J, "j"),
       col = cols, size = dat$value_adj, type = 's', aspect = c(1, 1, 1),
       xlab = "", ylab = "", zlab = "", axes = FALSE, box = TRUE)
aspect3d("iso")
vlab <- vlevels
dlab <- dlevels
jlab <- jlevels
axis3d("x-+", at = seq(vlab), labels = vlab, nticks = 2, cex.lab = 0.1)
axis3d("y-+", at = seq(dlab), labels = dlab, nticks = 2, cex.lab = 0.1)
axis3d("z+-", at = seq(jlab), labels = jlab, nticks = 2, cex.lab = 0.1)


grid3d("x--", at = seq(vlab))
grid3d("y--", at = seq(vlab))
grid3d("z-+", at = seq(vlab))

par3d(windowRect = c(75, 76, 1297, 749))
rgl.viewpoint(theta = 35, phi = 40, fov = 0, zoom = 0.7)

#rgl.postscript(paste0(fig.name, ".3d.pdf"), "pdf")
rgl.snapshot(paste0(fig.name, ".3d.png"))


#png(paste0(fig.name, ".color.png"), 230, 830, res = 100)
cols.std <- .colramp(0:15, prob)
pdf(paste0(fig.name, ".color.pdf"), 2, 8)
par(mar = c(2, 2, 2, 2))
legend_cols <- as.raster(matrix(rev(cols.std), ncol = 1))
plot(c(0, 2), c(0, 15), type = "n", axes = FALSE, xlab = "", ylab = "", main = "")
if (prob == "TRUE") {
    text(x = 1.6, y = c(0, 2, 4, 6, 8, 10, 15),
         labels = paste0(c(0, 2, 4, 6, 8, 10, 15), "%"), offset = 0, adj = 1)
} else {
    text(x = 1.8, y = seq(0, 1, l = 5), labels = round(seq(0, 1, l = 5) * max(dat$value)), offset = 0, adj = 1)
}
rasterImage(legend_cols, 0, 0, 1, 15)
dev.off()







