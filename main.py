import matplotlib.pyplot as plt

# Global constants ~ US
total_pop = 331449281  # Total population
beta = 0.3  # Infection rate
nu = 0.1  # Recovery Rate
initial_infected = 7
total_pop -= initial_infected


def dS_dt(S, I, beta):
  return -beta * I * S / total_pop


def dI_dt(I, S, beta, nu):
  return beta * I * S / total_pop - nu * I


def dR_dt(I, nu):
  return nu * I


def eulers_method(h, n, beta, nu):
  S = [total_pop]  # Initial susceptible population
  I = [initial_infected]  # Initial infected population
  R = [0]  # Initial recovered population
  t = [0]  # Time steps

  for i in range(1, n + 1):
    S.append(S[i - 1] + h * dS_dt(S[i - 1], I[i - 1], beta))
    I.append(I[i - 1] + h * dI_dt(I[i - 1], S[i - 1], beta, nu))
    R.append(R[i - 1] + h * dR_dt(I[i - 1], nu))
    t.append(t[i - 1] + h)

  return t, S, I, R


# Parameters
h = 1 / 365  # Step size
n = int((3 * 365) / h)  # Number of steps

# Running the Euler's method
t, S, I, R = eulers_method(h, n, beta, nu)

# Plotting
plt.plot(t, S, label='Susceptible')
plt.plot(t, I, label='Infectious')
plt.plot(t, R, label='Recovered')
plt.xlabel('Time (days)')
plt.ylabel('Population (hundered million)')
plt.title('USA: SIR Model')
plt.legend()
plt.show()
