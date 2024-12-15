from flask import Flask, render_template, request, jsonify
import math
import json

app = Flask(__name__)

# Load dataset
with open('constants.json') as f:
    constants = json.load(f)

@app.route("/")
def index():
    return render_template("index.html", constants=constants)

# Generate suggestions for autocomplete
suggestions = []
for category in constants.values():
    if isinstance(category, dict):
        suggestions.extend(category.keys())

@app.route('/autocomplete', methods=['GET'])
def autocomplete():
    query = request.args.get('q', '').lower()
    print(f"Query received: {query}")  # Debugging line

    if query:
        matches = [s for s in suggestions if query in s.lower()]
    else:
        matches = suggestions

    print(f"Matches found: {matches}")  # Debugging line
    return jsonify(matches)

@app.route("/calculate", methods=["POST"])
def calculate():
    data = request.json
    concentration = float(data.get("concentration", 0))  # Molar concentration
    compound = data.get("compound", "").lower()
    volume = float(data.get("volume", 1))  # Volume in Liters (default is 1L)

    if concentration <= 0:
        return jsonify({"error": "Concentration must be greater than 0."}), 400

    # Polyprotic Acids
    if "polyprotic_acids" in constants and compound in constants["polyprotic_acids"]:
        Ka_values = constants["polyprotic_acids"][compound]
        Ka_keys = sorted([key for key in Ka_values.keys() if key.startswith('Ka')], key=lambda k: int(k[2:]))

        H_plus_total = 0
        previous_H_plus = 0
        for key in Ka_keys:
            Ka = Ka_values[key]
            # Calculate [H+] for this dissociation step
            current_H_plus = math.sqrt(max(Ka * (concentration + previous_H_plus), 1e-14))
            H_plus_total += current_H_plus
            previous_H_plus = current_H_plus

        pH = -math.log10(max(H_plus_total, 1e-14))
        return jsonify({
            "pH": round(pH, 2),
            "pOH": round(14 - pH, 2)
        })

    # Polyprotic Bases
    if "polyprotic_bases" in constants and compound in constants["polyprotic_bases"]:
        Kb_values = constants["polyprotic_bases"][compound]
        Kb_keys = sorted([key for key in Kb_values.keys() if key.startswith('Kb')], key=lambda k: int(k[2:]))

        OH_minus_total = 0
        previous_OH_minus = 0
        for key in Kb_keys:
            Kb = Kb_values[key]
            # Calculate [OH-] for this dissociation step
            current_OH_minus = math.sqrt(max(Kb * (concentration + previous_OH_minus), 1e-14))
            OH_minus_total += current_OH_minus
            previous_OH_minus = current_OH_minus

        pOH = -math.log10(max(OH_minus_total, 1e-14))
        return jsonify({
            "pH": round(14 - pOH, 2),
            "pOH": round(pOH, 2)
        })

    # Strong Acids and Bases
    for category, compounds in constants.items():
        if compound in compounds:
            dissociation_constant = compounds[compound]

            # Handle strong acids
            if category == "strong_acids":
                pH = -math.log10(max(concentration, 1e-14))  # Avoid log10(0)
                return jsonify({
                    "pH": round(pH, 2),
                    "pOH": round(14 - pH, 2)
                })

            # Handle strong bases
            if category == "strong_bases":
                pOH = -math.log10(max(concentration, 1e-14))  # Avoid log10(0)
                return jsonify({
                    "pH": round(14 - pOH, 2),
                    "pOH": round(pOH, 2)
                })

            # Handle weak acids
            if category == "weak_acids":
                H_plus = math.sqrt(max(dissociation_constant * concentration, 1e-14))
                pH = -math.log10(max(H_plus, 1e-14))
                return jsonify({
                    "pH": round(pH, 2),
                    "pOH": round(14 - pH, 2)
                })

            # Handle weak bases
            if category == "weak_bases":
                OH_minus = math.sqrt(max(dissociation_constant * concentration, 1e-14))
                pOH = -math.log10(max(OH_minus, 1e-14))
                return jsonify({
                    "pH": round(14 - pOH, 2),
                    "pOH": round(pOH, 2)
                })

    return jsonify({"error": "Compound not found."}), 400

if __name__ == "__main__":
    app.run(debug=True)
