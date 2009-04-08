library(qtl)
f2qtl <- c(6,150,3,7)
map <- sim.map(len=rep(200,10), n.mar=10, eq.spacing=T)
f2cross <- sim.cross(map,f2qtl,n=100,type="f2")       # Simulate a F2 Cross
cof <- MQMCofactorsEach(f2cross,8)
aa <- scanMQM(f2cross,cof)
