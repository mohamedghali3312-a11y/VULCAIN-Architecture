# VULCAIN (v5.1) — Architecture Hybride d'Algèbre Exacte en Mémoire O(1) et de Factorisation Polynomiale Certifiée

**Statut :** `CI : PASSING` | **HAL :** `hal-05860758` | **DOI :** `10.5281/zenodo.24311057` | **Preuves Formelles :** `Lean 4 (Mathlib)` | **Licence :** `MIT (2026)`

**Auteur :** Mohamed Ghali — Chercheur indépendant (France) | Novembre 2026

---

## 1. Vue d'ensemble de la Trilogie VULCAIN

**VULCAIN (v5.1)** unifie en un pipeline hybride à trois étages l'algèbre modulaire exacte en mémoire constante et l'analyse numérique certifiée à haute précision :

* **Étage 1 — L'Enclume Discrète (*La Machine v16.0*) :** Table circulaire miroir V₄ = 2 310 en cache CPU L1 (< 2 Ko), filtre PGCD W₅..₁₂ (85,49 % d'élimination), transformée **NTT-Montgomery** sans division, Théorème des Restes Chinois (CRT) et certification déterministe 2-adique prouvée sous **Lean 4** (`DOI: 10.5281/zenodo.24334996`).
* **Étage 2 — Le Recuit Continu (*MBC-Annealing v5*) :** Factorisation simultanée en blocs quadratiques réels Bⱼ(x) = x² + Sⱼx + Pⱼ régularisés sur le cône semi-défini positif avec arithmétique compensée **EFT (106 bits)** et certification a posteriori de **Smale-Kantorovich** (`DOI: 10.5281/zenodo.24284263`).
* **Étage 3 — Le Raffinement Vectorisé (*VULCAIN v5.1*) :** Passe finale vectorisée d'**Aberth-Ehrlich** (*overflow-safe*) atteignant la précision machine IEEE 754 (`DOI: 10.5281/zenodo.24311057`).

---

## 2. Citation académique officielle

> Ghali, M. (2026). *VULCAIN (v5.1): Hybrid Architecture for O(1)-Memory Exact Algebra and Certified High-Precision Polynomial Factorization*. HAL (`hal-05860758`) / Zenodo (`DOI: 10.5281/zenodo.24311057`).

---

## 3. Structure du dépôt

```text
VULCAIN-Architecture/
├── .github/workflows/ci.yml      # Intégration continue automatique (Python + Lean 4)
├── vulcain_v5_complete.py        # Moteur hybride complet à 3 étages (NTT-CRT + EFT + Aberth)
├── test_vulcain.py               # Suite de tests unitaires (pytest)
├── proofs_vulcain.lean           # Preuves formelles certifiées sous Lean 4
├── lakefile.toml                 # Configuration Lean 4 / Mathlib
├── lake-manifest.json            # Verrouillage des dépendances Lean 4
├── lean-toolchain                # Version du compilateur Lean (v4.12.0)
└── LICENSE                       # Licence MIT (Mohamed Ghali, 2026)
```

---

## 4. Exécution et vérification locale

```bash
# Tests unitaires et démonstration du pipeline hybride Python
pip install pytest numpy
pytest -v
python vulcain_v5_complete.py

# Vérification formelle des théorèmes sous Lean 4
lake update
lake exe cache get
lake build
```
