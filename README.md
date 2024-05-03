# Edward Summation

## NOTE

1. Con σ → ∞, si torna al caso classico. Le due distribuzioni di densità gaussiane tendono a una distribuzione costante di intensità nulla. Sotto queste condizioni e l'energia totale è conservata. Il sistema si comporta come il caso classico. Con σ infinita conta solo la componente spazio reale, con sigma nulla conta solo la componente di spazio reciproco. Perchè con σ nulla le particelle nello spazio reale sono equivalenti a due δ(0) sovrapposte, una con -Q e una con +Q, la carica risultante è nulla.

2. Per σ ≠ ∞ l'energia non è conservata. Dovrebbe esserlo? Possibile errore nell'implementazione della componente spazio reciproco della edwald summation.

3. Test effettuato per 270 particelle e ALPHA = 5.6/CELL_L. Risultati (Paper 2 sezione 3.3):

   1. REAL SUM: ≈1.2e-3
   2. RECIPROCAL SUM: ≈-1.9e-4

4. La somma delle forze reciproche è negativa perchè è dovuta a gausiane della stessa carica della particella, tende a schiacciare la particella verso l'origine.

5. La componente di spazio reciproco non dipende da r, anche se la carica non è nulla nella scatola (0,0,0) posso ignorare il suo contributo.
