"""
MONTE CARLO OPTION PRICING
===========================
Pricing di opzioni finanziarie tramite simulazione Monte Carlo di un moto
browniano geometrico (GBM), sotto la misura di probabilità risk-neutral.

1. Call europea: payoff dipende solo dal prezzo finale S_T.
   Esiste formula chiusa (Black-Scholes), usata per VALIDARE la simulazione.
2. Opzione asiatica: payoff dipende dalla media del prezzo sull'intera
   traiettoria (path-dependent). Non esiste formula chiusa: qui Monte
   Carlo e' lo strumento di calcolo standard, non solo una verifica.
"""

import numpy as np
from scipy.stats import norm

np.random.seed(42)

# --- Parametri comuni ---
S_0 = 100     # prezzo iniziale
K = 105       # strike price
r = 0.03      # tasso privo di rischio
sigma = 0.20  # volatilita' annua
T = 1         # scadenza, in anni

n_simulazioni = 100_000
n_passi = 252
dt = T / n_passi
time_grid = np.linspace(0, T, n_passi + 1)  # vettore dei tempi


# ============================================================
# PARTE 1 - Call europea (payoff dipende solo da S_T)
# ============================================================

Z = np.random.normal(loc=0, scale=1, size=n_simulazioni)

S_T = S_0 * np.exp((r - 0.5 * sigma**2) * T + sigma * np.sqrt(T) * Z)

payoff_call = np.maximum(S_T - K, 0)
V_0_call_mc = np.exp(-r * T) * payoff_call.mean()

# Formula chiusa di Black-Scholes, per validare la simulazione
d1 = (np.log(S_0 / K) + (r + 0.5 * sigma**2) * T) / (sigma * np.sqrt(T))
d2 = d1 - sigma * np.sqrt(T)
V_0_call_bs = S_0 * norm.cdf(d1) - K * np.exp(-r * T) * norm.cdf(d2)

print("=== Call europea ===")
print(f"Prezzo Monte Carlo:   {V_0_call_mc:.4f}")
print(f"Prezzo Black-Scholes: {V_0_call_bs:.4f}")
print(f"Errore relativo:      {abs(V_0_call_mc - V_0_call_bs) / V_0_call_bs * 100:.3f}%")


# ============================================================
# PARTE 2 - Opzione asiatica (payoff dipende dalla media sulla traiettoria)
# ============================================================

incrementi_browniani = np.random.normal(0, np.sqrt(dt), size=(n_simulazioni, n_passi))
brownian_paths = np.concatenate(
    [np.zeros((n_simulazioni, 1)), np.cumsum(incrementi_browniani, axis=1)], axis=1
)
# matrice con tutte le traiettorie browniane

S_t = S_0 * np.exp((r - 0.5 * sigma**2) * time_grid + sigma * brownian_paths)
# matrice con tutte le traiettorie dei prezzi
# righe = simulazioni/traiettorie diverse
# colonne = tempi

prezzi_medi = np.mean(S_t, axis=1)
payoff_asiatica = np.maximum(prezzi_medi - K, 0)

expected_montecarlo_payoff = np.mean(payoff_asiatica)
V_0_asiatica = np.exp(-r * T) * expected_montecarlo_payoff

print("\n=== Opzione asiatica (average price) ===")
print(f"Prezzo Monte Carlo:   {V_0_asiatica:.4f}")
print("(nessuna formula chiusa semplice disponibile per confronto)")

print(f"\nL'opzione asiatica costa meno della call europea: {V_0_asiatica < V_0_call_mc} "
      f"(atteso: mediare sulla traiettoria riduce la varianza del payoff)")
