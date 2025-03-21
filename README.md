# Sequence Analyzer

🔬 **Sequence Analyzer** is a bioinformatics tool designed to simplify sequence analysis and alignment tasks. Built using **Streamlit**, **Biopython**, and **Plotly**, this application provides an intuitive interface for analyzing DNA and RNA sequences, aligning sequences, and fetching sequences from GenBank.

---

## 🚀 Features

### 1. **Sequence Alignment**

- Perform **Pairwise Alignment** or **Multiple Sequence Alignment (MSA)** for DNA or RNA sequences.
- Upload sequences in FASTA format or fetch them directly from GenBank using accession numbers.
- Download aligned sequences in `.fasta`, `.txt`, or `.rtf` formats.

### 2. **DNA Analysis**

- Analyze DNA sequences for:
  - GC content
  - Nucleotide composition (A, T, G, C)
- Visualize results with interactive histograms and scatter plots.
- Export analysis results in `.csv`, `.txt`, or `.rtf` formats.

### 3. **RNA Analysis**

- Analyze RNA sequences for:
  - GC content
  - Nucleotide composition (A, U, G, C)
- Visualize results with interactive histograms and scatter plots.
- Export analysis results in `.csv`, `.txt`, or `.rtf` formats.

### 4. **GenBank Sequence Retrieval**

- Fetch sequences from GenBank using accession numbers.
- Supports integration with NCBI's Entrez API.

### 5. **Contact and About Pages**

- Learn more about the application and its contributors.
- Contact information for inquiries and contributions.

---

## 📂 Project Structure

```txt
sequence_analyzer/
│
├── app.py                     # Main application entry point
├── utils.py                   # Utility functions for sequence analysis and visualization
├── requirements.txt           # Python dependencies
├── pages/                     # Streamlit multipage application
│   ├── 1_Sequence_Alignment.py # Sequence alignment page
│   ├── 2_DNA_Analysis.py       # DNA analysis page
│   ├── 3_RNA_Analysis.py       # RNA analysis page
│   ├── 4_About.py              # About page
│   └── 5_Contact.py            # Contact page
├── README.md                  # Project documentation
└── .gitignore                 # Git ignore file
```

---

## 🛠️ Installation

1. **Clone the Repository**

   ```bash
   git clone https://github.com/bioinformatics-project/sequence_analyzer.git
   cd sequence_analyzer
   ```
2. **Set Up a Virtual Environment**

   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```
3. **Install Dependencies**

   ```bash
   pip install -r requirements.txt
   ```
4. **Run the Application**

   ```bash
   streamlit run app.py
   ```

---

## 📊 Usage

1. Launch the application using the command above.
2. Use the **navigation menu** on the left to explore the following features:
   - **Sequence Alignment**: Upload or fetch sequences and perform alignments.
   - **DNA Analysis**: Analyze DNA sequences for GC content and nucleotide composition.
   - **RNA Analysis**: Analyze RNA sequences for GC content and nucleotide composition.
   - **About**: Learn more about the application.
   - **Contact**: Reach out for inquiries or contributions.

---

## 📜 License

This project is licensed under the [MIT License](LICENSE).

---

## 🤝 Contributing

We welcome contributions to improve this application! To contribute:

1. Fork the repository.
2. Create a new branch for your feature or bug fix.
3. Submit a pull request with a detailed description of your changes.

---

## 📧 Contact

For inquiries or contributions, please reach out:

- 📧 Email: **[Muhammad Abiodun SULAIMAN](mailto:abiodun.msulaiman@gmail.com), [Bolaji Fatai OYEYEMI](mailto:bolajioyeyemi@gmail.com)**
- 🛠 GitHub: **[Sequence Analyzer](https://github.com/Behordeun/sequence_analyzer)**

---

## ❤️ Acknowledgments

Developed with ❤️ using **Streamlit**, **Biopython**, and **Plotly** by:

- [Behordeun](https://github.com/Behordeun)
- [Bollergene](https://github.com/bollergene)
