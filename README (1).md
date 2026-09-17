# Monte Carlo Option Pricing

Pricing di opzioni finanziarie tramite simulazione Monte Carlo di un moto browniano geometrico (GBM), con validazione contro la formula chiusa di Black-Scholes e applicazione a un payoff path-dependent per cui non esiste soluzione analitica in forma chiusa.

## Motivazione

Il metodo Monte Carlo stima un valore atteso `E[X]` generando `N` realizzazioni indipendenti e identicamente distribuite `X_1, ..., X_N` della variabile aleatoria `X`, e approssimandolo con la media campionaria:

```
E[X] ≈ (1/N) Σ X_i
```

La validità di questa approssimazione, per `N → ∞`, è garantita dalla **Legge dei Grandi Numeri (LLN)**.

In finanza quantitativa, sotto la misura di probabilità **risk-neutral** `Q` — il cui utilizzo è giustificato dal **teorema di Girsanov** (cambio di misura assolutamente continuo che elimina il drift di mercato a favore del tasso privo di rischio) — il prezzo `V_0` di un derivato con payoff `Φ` a scadenza `T` è il valore atteso, sotto `Q`, del payoff scontato:

```
V_0 = exp(-r·T) · E^Q[ Φ ]
```

Questo progetto nasce dall'applicazione pratica di concetti approfonditi nella tesi di laurea magistrale in Matematica ("The Girsanov Theorem with Applications to the Theory of SDEs"), traducendo in codice la teoria delle equazioni differenziali stocastiche (SDE) e del cambio di misura risk-neutral.

Si considerano due casi:

1. **Call europea** — payoff `Φ = max(S_T - K, 0)`, dipendente solo dal prezzo a scadenza. Esiste una formula chiusa (Black-Scholes): la si usa per **validare** la simulazione Monte Carlo.
2. **Opzione asiatica (average price)** — payoff dipendente dalla media del prezzo sull'intera traiettoria simulata, non solo dal valore a scadenza. Per questo payoff path-dependent non esiste una formula chiusa semplice: qui Monte Carlo non è un'alternativa, è lo strumento di calcolo standard usato nella pratica.

## Notazione

| Simbolo | Significato |
|---|---|
| `S_0` | Prezzo del sottostante al tempo 0 |
| `S_t` | Prezzo del sottostante al tempo `t` |
| `K` | Strike price dell'opzione |
| `r` | Tasso di interesse privo di rischio (risk-free rate) |
| `σ` | Volatilità (deviazione standard istantanea dei log-rendimenti) |
| `T` | Scadenza dell'opzione (maturity), in anni |
| `W_t` | Moto browniano standard (processo di Wiener) sotto la misura `Q` |
| `N` | Numero di traiettorie simulate |

## Metodologia

Sotto la misura risk-neutral `Q`, il prezzo del sottostante segue la SDE:

```
dS_t = r·S_t·dt + σ·S_t·dW_t
```

Applicando il Lemma di Itô a `Y_t = ln(S_t)`, si ottiene la soluzione esplicita del GBM:

```
S_t = S_0 · exp( (r - σ²/2)·t + σ·W_t )
```

Per ciascuna delle `N` traiettorie simulate si genera `W_t` come somma cumulata di incrementi gaussiani indipendenti `ΔW ~ N(0, Δt)`, si ottiene `S_t` applicando la formula sopra istante per istante, e si calcola il payoff di interesse sulla traiettoria completa.

## Risultati

Parametri: `S_0 = 100`, `K = 105`, `r = 3%`, `σ = 20%`, `T = 1` anno, `N = 100.000` traiettorie simulate.

| Opzione | Prezzo Monte Carlo | Prezzo teorico (Black-Scholes) | Errore relativo |
|---|---|---|---|
| Call europea | 7.1478 | 7.1281 | 0.28% |
| Asiatica (average price) | 3.1250 | — (nessuna formula chiusa) | — |

L'opzione asiatica risulta più economica della call europea, coerentemente con la teoria: mediare il prezzo sull'intera traiettoria riduce la varianza del payoff rispetto a considerare solo il valore a scadenza.

## Come eseguirlo

```bash
pip install numpy scipy
python monte_carlo_pricing.py
```

## Struttura del codice

Lo script e' organizzato in due parti sequenziali:

- **Parte 1 - Call europea**: simulazione diretta di `S_T` (formula esplicita del GBM applicata al tempo finale `T`), pricing Monte Carlo e confronto con la formula chiusa di Black-Scholes.
- **Parte 2 - Opzione asiatica**: costruzione della traiettoria completa del moto browniano (somma cumulata di incrementi gaussiani su `n_passi` intervalli di ampiezza `dt`), applicazione della formula del GBM istante per istante, media della traiettoria e pricing del payoff path-dependent.
