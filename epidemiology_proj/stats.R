#statistical tests

mse <- function(model, datum) {
  var = 0
  t = length(model)
  for (i in 1: t){
    var = var + (model[i]- datum[i])^2
  }
  return(var/t)
}

find_start <- function(datum, tmin, tmax) {
  for (t in tmin:tmax){
    s2 = datum$New.Cases[(t+1):(t+5)]
    s1 = datum$New.Cases[(t-5):(t-1)]
    x = t.test(s2, s1, alternative = "greater")
    if(x$p.value < 0.02){
      return(t)
    }
  }
  return(0)
}

num_susceptibles_sir <- function (datum, tmin, tmax, beta, gamma, psi1, psi2, psi3, nu){
  #determine I0, R0, and S0 at time tmin-1 that carries over into day tmin
  I0 = sum(datum$New.Cases[(tmin-9):(tmin)])
  P0 = nu * (datum$Total.Cases[(tmin-10)] - Total.Deaths[(tmin-10)])
  
  #Vaccines assumed to be effective in 2 weeks
  V0 = psi1*(datum$People.Receiving.1.or.More.Doses.Cumulative[(tmin-14)] -
               datum$People.Fully.Vaccinated.Cumulative[(tmin-14)]) +
    psi2*(datum$People.Fully.Vaccinated.Cumulative[(tmin-14)] -
            datum$People.Receiving.a.Booster.Dose.Cumulative[(tmin-14)]) + 
    psi3 * datum$People.Receiving.a.Booster.Dose.Cumulative[(tmin-14)]
  
  R0 = P0 + V0 *(1 - P0/N)
  S0 = N - I0 - R0 
  t_len = tmax - tmin
  
  S = vector()
  R = vector()
  I = vector()
  P = vector()
  
  dose1 = datum$People.Receiving.1.or.More.Doses.Cumulative[(tmin - 13):(tmax - 13)]
  dose2 = datum$People.Fully.Vaccinated.Cumulative[(tmin - 13):(tmax -13)]
  dose3 = datum$People.Receiving.a.Booster.Dose.Cumulative[(tmin - 13):(tmax -13)]
  
  V = psi1* (dose1 - dose2) + psi2* (dose2 - dose3) + psi3 * dose3
  S[1] = S0
  R[1] = R0
  I[1] = I0
  P[1] = P0
  
  #returned vector
  cumulative_cases = vector()
  cumulative_cases[1] = datum$Total.Cases[tmin]
  
  for (t in 1:t_len) {
    dI = beta*I[t]*S[t]/N - gamma*I[t]
    cumulative_cases[t+1] = cumulative_cases[t] + beta*I[t]*S[t]/N
    dP = gamma*I[t]
    prop = P[t]/N
    I[t+1]= I[t] + dI
    P[t+1] = P[t] + dP
    R[t+1]= P[t] + dP + V[t] *(1-prop)
    S[t+1]= N - I[t+1] - R[t+1]
    
  }
  return (c(S0, S[(t_len + 1)]))
}

num_susceptibles_sir_omi<- function (datum, tmin, tmax, beta, gamma, psi1, psi2, psi3, nu){
  I0 = sum(New.Cases[(tmin-9):(tmin)])
  P0 = nu * (Total.Cases[(tmin-10)] - Total.Deaths[(tmin-10)])
  V0 = psi1*(datum$People.Fully.Vaccinated.Cumulative[(tmin-28)] - 
               datum$People.Receiving.a.Booster.Dose.Cumulative[(tmin-14)]) +
    psi2*(datum$People.Fully.Vaccinated.Cumulative[(tmin-14)] -
            datum$People.Fully.Vaccinated.Cumulative[(tmin-28)]) + 
    psi3 * datum$People.Receiving.a.Booster.Dose.Cumulative[(tmin-14)]
  
  R0 = P0 + V0*(1 - P0/N)
  S0 = N - I0 - R0 
  t_len = tmax - tmin
  
  S = vector()
  R = vector()
  I = vector()
  P = vector()
  
  dose2_old = datum$People.Fully.Vaccinated.Cumulative[(tmin - 27):(tmax - 27)]- 
    datum$People.Receiving.a.Booster.Dose.Cumulative[(tmin-13):(tmax - 13)]
  dose2_new = datum$People.Fully.Vaccinated.Cumulative[(tmin - 13):(tmax -13)] - 
    datum$People.Fully.Vaccinated.Cumulative[(tmin-27):(tmax-27)]
  dose3 = datum$People.Receiving.a.Booster.Dose.Cumulative[(tmin - 13):(tmax -13)]
  
  V = psi1* dose2_old + psi2* dose2_new + psi3 * dose3
  
  S[1] = S0
  R[1] = R0
  I[1] = I0
  P[1] = P0
  
  #returned vector
  cumulative_cases = vector()
  cumulative_cases[1] = datum$Total.Cases[tmin]
  
  for (t in 1:t_len) {
    dI = beta*I[t]*S[t]/N - gamma*I[t]
    cumulative_cases[t+1] = cumulative_cases[t] + beta*I[t]*S[t]/N
    dP = gamma*I[t]
    prop = P[t]/N
    I[t+1]= I[t] + dI
    P[t+1] = P[t] + dP
    R[t+1]= P[t] + dP + V[t]*(1-prop)
    S[t+1]= N - I[t+1] - R[t+1]
  }
  return(c(S0, S[(t_len + 1)]))
}
