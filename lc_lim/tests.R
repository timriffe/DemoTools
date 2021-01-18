

input <- bind_rows(
  readxl::read_xlsx("lc/Li_2018_Limited_Lee-Carter-v4.xlsm",sheet = "Input_M",range = "A1:F23") %>% 
    rename(Age = Input_M) %>% mutate(Sex = "m"),
  readxl::read_xlsx("lc/Li_2018_Limited_Lee-Carter-v4.xlsm",sheet = "Input_F",range = "A1:F23") %>% 
    rename(Age = Input_F) %>% mutate(Sex = "f")) %>% 
  select(Sex, Age, everything())
dput(as.matrix(round(input[23:44,-c(1:2)],6)))


# troubles ----------------------------------------------------------------

Age <- c(0,1,seq(5,100,5))
Years <- as.double(sort(colnames(input[-c(1:2)])))
Years_target <- 1948 + 1:14*5
nYears_Target <- length(Years_target)
e0_Males <- 75:78
e0_Females <- 77:80
e0_Years <- seq(1980,2010,10)
ne0_Years = length(e0_Years)
prev_divergence = F
Males <- as.matrix(input[1:22,-c(1:2)])
Females <- as.matrix(input[23:44,-c(1:2)])
type = "m"
extrapLaw = "Coale-Guo"


library(DemoTools)
age = c(0,1,seq(5,100,5))
m1 = c(0.105,0.0063,7e-04,6e-04,9e-04,0.0011,0.0015,0.002,0.003,0.0041,0.0048,0.0077,0.0098,0.0156,0.025,0.0317,0.0572,0.0597,0.0663,0.0835,0.1192,0.193)
m2 = c(0.0824,0.0048,7e-04,5e-04,8e-04,0.001,0.0014,0.0019,0.0028,0.0038,0.0045,0.0072,0.0092,0.0146,0.0231,0.03,0.0527,0.0575,0.0664,0.0859,0.1245,0.2021)
m3 = c(0.0647,0.0037,6e-04,5e-04,8e-04,0.001,0.0013,0.0017,0.0025,0.0034,0.0042,0.0066,0.0087,0.0137,0.0214,0.0284,0.0485,0.0554,0.0665,0.0884,0.13,0.2116)
plot(age,m1,log = "y"); lines(age,m2,log = "y",col=2); lines(age,m3,log = "y",col=3)
plot( lt_abridged(nMx = m1, Age = age, Sex = "f")$ex, ylim=c(0,80))
lines(lt_abridged(nMx = m2, Age = age, Sex = "f")$ex,col=2)
lines(lt_abridged(nMx = m3, Age = age, Sex = "f")$ex,col=3)


m1 = nMxf_hat[,3]
m2 = nMxf_hat[,4]
m3 = nMxf_hat[,5]

library(DemoTools)
age = c(0,1,seq(5,100,5))
m1 = c(0.105,0.0063,7e-04,6e-04,9e-04,0.0011,0.0015,0.002,0.003,0.0041,0.0048,0.0077,0.0098,0.0156,0.025,0.0317,0.0572,0.0597,0.0663,0.0835,0.1192,0.193)
m2 = c(0.0824,0.0048,7e-04,5e-04,8e-04,0.001,0.0014,0.0019,0.0028,0.0038,0.0045,0.0072,0.0092,0.0146,0.0231,0.03,0.0527,0.0575,0.0664,0.0859,0.1245,0.2021)
m3 = c(0.0647,0.0037,6e-04,5e-04,8e-04,0.001,0.0013,0.0017,0.0025,0.0034,0.0042,0.0066,0.0087,0.0137,0.0214,0.0284,0.0485,0.0554,0.0665,0.0884,0.13,0.2116)
plot(age,m1,log = "y"); lines(age,m2,log = "y",col=2); lines(age,m3,log = "y",col=3)
plot( lt_abridged(nMx = m1, Age = age, Sex = "f")$ex, ylim=c(0,80))
lines(lt_abridged(nMx = m2, Age = age, Sex = "f")$ex,col=2)
lines(lt_abridged(nMx = m3, Age = age, Sex = "f")$ex,col=3)


plot( lt_abridged(nMx = m1, Age = age, Sex = "f", axmethod = "un")$ex, ylim=c(0,80))
lines(lt_abridged(nMx = m2, Age = age, Sex = "f", axmethod = "un")$ex,col=2)
lines(lt_abridged(nMx = m3, Age = age, Sex = "f", axmethod = "un")$ex,col=3)

plot( lt_abridged(nMx = m1, Age = age, Sex = "f", a0rule = "cd")$ex, ylim=c(0,80))
lines(lt_abridged(nMx = m2, Age = age, Sex = "f", a0rule = "cd")$ex,col=2)
lines(lt_abridged(nMx = m3, Age = age, Sex = "f", a0rule = "cd")$ex,col=3)



#######################

nMxm = Males

axm = rowSums(log(nMxm))/nYears
ktom = colSums(log(nMxm))-sum(axm)
bxm = rowSums(sweep(log(nMxm) - axm, MARGIN = 2, ktom, `*`))/sum(ktom^2)

lc_lim_data %>% 
  filter(Age==0) %>% select(Year,Sex,ex) %>% 
  ggplot() + geom_line(aes(Year,ex,col=Sex)) +
  geom_line()










# interp/extrap
# interp e0
e0pm = c(e0m[1] + (ywpp[ywpp<ye0[1]] - ye0[1]) / (ye0[2] - ye0[1]) * (e0m[2] - e0m[1]),
         approxfun(x=ye0, y = e0m, method = "linear")(ywpp[ywpp>ye0[1] & ywpp<ye0[nye0]]),
         e0m[nye0] + (ywpp[ywpp>ye0[nye0]] - ye0[nye0]) / (ye0[nye0] - ye0[nye0-1]) * (e0m[nye0] - e0m[nye0-1]))
e0pf = c(e0f[1] + (ywpp[ywpp<ye0[1]] - ye0[1]) / (ye0[2] - ye0[1]) * (e0f[2] - e0f[1]),
         approxfun(x=ye0, y = e0f, method = "linear")(ywpp[ywpp>ye0[1] & ywpp<ye0[nye0]]),
         e0f[nye0] + (ywpp[ywpp>ye0[nye0]] - ye0[nye0]) / (ye0[nye0] - ye0[nye0-1]) * (e0f[nye0] - e0f[nye0-1]))
plot(ywpp,e0pm);lines(ye0,e0m,col=2)



