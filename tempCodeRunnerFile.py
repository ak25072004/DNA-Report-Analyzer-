from flask import Flask, render_template, request
import os

app = Flask(__name__)

# Utility functions
def load_user_dna(file_stream):
    return file_stream.read().decode("utf-8").strip()

def load_disease_markers(file_path):
    markers = {}
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 2:
                disease, marker = parts
                markers[disease] = marker
    return markers

def match_score(seq1, seq2):
    return sum(1 for a, b in zip(seq1, seq2) if a == b)

# Rabin-Karp algorithm
def rabin_karp_match(user_dna, marker, threshold=0.8):
    base = 256
    mod = 101
    m = len(marker)
    if len(user_dna) < m:
        return False

    h_marker = 0
    h_segment = 0
    h = 1

    for i in range(m - 1):
        h = (h * base) % mod

    for i in range(m):
        h_marker = (base * h_marker + ord(marker[i])) % mod
        h_segment = (base * h_segment + ord(user_dna[i])) % mod

    for i in range(len(user_dna) - m + 1):
        if h_marker == h_segment:
            if user_dna[i:i + m] == marker:
                return True
        if i < len(user_dna) - m:
            h_segment = (base * (h_segment - ord(user_dna[i]) * h) + ord(user_dna[i + m])) % mod
            h_segment = (h_segment + mod) % mod
    return False

# Longest Common Subsequence
def lcs(a, b):
    n, m = len(a), len(b)
    dp = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(n):
        for j in range(m):
            if a[i] == b[j]:
                dp[i + 1][j + 1] = dp[i][j] + 1
            else:
                dp[i + 1][j + 1] = max(dp[i][j + 1], dp[i + 1][j])
    return dp[n][m]

# Jaccard Similarity
def jaccard_similarity(a, b, k=5):
    a_kmers = set(a[i:i + k] for i in range(len(a) - k + 1))
    b_kmers = set(b[i:i + k] for i in range(len(b) - k + 1))
    intersection = a_kmers & b_kmers
    union = a_kmers | b_kmers
    if not union:
        return 0
    return len(intersection) / len(union)

# DNA disease comparison
def compare_dna(user_dna, markers, method="simple"):
    detected = []

    for disease, marker in markers.items():
        if method == "rabin":
            if rabin_karp_match(user_dna, marker):
                detected.append(disease)

        elif method == "lcs":
            max_score = max(
                lcs(user_dna[i:i + len(marker)], marker)
                for i in range(len(user_dna) - len(marker) + 1)
            )
            if max_score >= len(marker) * 0.8:
                detected.append(disease)

        elif method == "jaccard":
            max_score = max(
                jaccard_similarity(user_dna[i:i + len(marker)], marker)
                for i in range(len(user_dna) - len(marker) + 1)
            )
            if max_score >= 0.5:
                detected.append(disease)

        else:  # Simple match
            for i in range(len(user_dna) - len(marker) + 1):
                segment = user_dna[i:i + len(marker)]
                score = match_score(segment, marker)
                if score >= len(marker) * 0.8:
                    detected.append(disease)
                    break

    return detected

# DNA similarity comparison between two samples
def calculate_similarity(dna1, dna2, method="simple"):
    min_len = min(len(dna1), len(dna2))
    if method == "rabin":
        k = 10
        matches = 0
        total = min_len - k + 1
        for i in range(total):
            segment = dna1[i:i + k]
            if rabin_karp_match(dna2, segment):
                matches += 1
        return matches / total if total > 0 else 0

    elif method == "lcs":
        return lcs(dna1, dna2) / min_len if min_len > 0 else 0

    elif method == "jaccard":
        return jaccard_similarity(dna1, dna2)

    else:  # simple
        matches = sum(1 for a, b in zip(dna1, dna2) if a == b)
        return matches / min_len if min_len > 0 else 0

# Flask Route
@app.route("/", methods=["GET", "POST"])
def index():
    result = {}
    selected_method = "simple"

    if request.method == "POST":
        selected_method = request.form.get("method", "simple")
        file1 = request.files.get("dna1")
        file2 = request.files.get("dna2")

        if not file1 or not file2:
            return render_template("index.html", results=None, selected_method=selected_method)

        user_dna1 = load_user_dna(file1)
        user_dna2 = load_user_dna(file2)

        disease_markers = load_disease_markers("diseases.txt")

        risks1 = compare_dna(user_dna1, disease_markers, method=selected_method)
        risks2 = compare_dna(user_dna2, disease_markers, method=selected_method)
        similarity_score = calculate_similarity(user_dna1, user_dna2, selected_method)
        similarity_percent = f"{similarity_score * 100:.2f}%"

        result = {
            "sample1": risks1,
            "sample2": risks2,
            "similarity": similarity_percent
        }

    return render_template("index.html", results=result, selected_method=selected_method)

if __name__ == "__main__":
    app.run(debug=True)