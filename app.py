from flask import Flask, render_template, request, jsonify
import math
import json

app = Flask(__name__)

# Load dataset
with open('constants.json') as f:
    constants = json.load(f)

# Extract possible compounds from constants for autocomplete suggestions
suggestions = []
for category in constants.values():
    if isinstance(category, dict):
        suggestions.extend(category.keys())

@app.route("/")
def index():
    return render_template("index.html")

@app.route('/autocomplete', methods=['GET'])
def autocomplete():
    query = request.args.get('q', '').lower()
    if query:
        matches = [s for s in suggestions if query in s.lower()]
        return jsonify(matches[:10])  # Limit to 10 suggestions
    return jsonify([])

@app.route("/calculate", methods=["POST"])
def calculate():
    try:
        compound = request.form.get("compound", "").lower()
        calculation_mode = request.form.get("calculation_mode", "")
        
        # Find compound data and category
        compound_data = None
        compound_category = None
        for category, compounds in constants.items():
            if compound in compounds:
                compound_data = compounds[compound]
                compound_category = category
                break
        
        if not compound_data:
            return jsonify({"error": "Compound not found"}), 400

        # Calculate concentration based on mode
        if calculation_mode == "concentration":
            concentration = float(request.form.get("concentration", 0))
            if concentration <= 0:
                return jsonify({"error": "Concentration must be greater than 0"}), 400

        elif calculation_mode in ["moles-volume", "mass-volume"]:
            volume = float(request.form.get("volume", 0))
            if volume <= 0:
                return jsonify({"error": "Volume must be greater than 0"}), 400

            if calculation_mode == "moles-volume":
                moles = float(request.form.get("moles", 0))
                if moles <= 0:
                    return jsonify({"error": "Moles must be greater than 0"}), 400
                concentration = moles / volume

            else:  # mass-volume mode
                mass = float(request.form.get("mass", 0))
                if mass <= 0:
                    return jsonify({"error": "Mass must be greater than 0"}), 400
                
                molecular_weight = compound_data.get("Mw", 0)
                if molecular_weight <= 0:
                    return jsonify({"error": "Molecular weight not available for this compound"}), 400
                
                moles = mass / molecular_weight
                concentration = moles / volume

        else:
            return jsonify({"error": "Invalid calculation mode"}), 400

        # Calculate pH and pOH based on compound type
        if compound_category == "strong_acids":
            pH = -math.log10(concentration)
            pOH = 14 - pH

        elif compound_category == "strong_bases":
            pOH = -math.log10(concentration)
            pH = 14 - pOH

        elif compound_category == "weak_acids":
            Ka = compound_data["K"]
            H_plus = math.sqrt(Ka * concentration)
            pH = -math.log10(H_plus)
            pOH = 14 - pH

        elif compound_category == "weak_bases":
            Kb = compound_data["K"]
            OH_minus = math.sqrt(Kb * concentration)
            pOH = -math.log10(OH_minus)
            pH = 14 - pOH

        elif compound_category == "polyprotic_acids":
            # Handle multi-step dissociation
            Ka_values = []
            for key in sorted(compound_data["K"].keys()):  # Sort to ensure correct order
                if key.startswith(('Ka', 'step_')):
                    Ka_values.append(compound_data["K"][key])
            
            H_plus_total = 0
            for Ka in Ka_values:
                H_plus = math.sqrt(Ka * concentration)
                H_plus_total += H_plus
            
            pH = -math.log10(H_plus_total)
            pOH = 14 - pH

        elif compound_category == "polyprotic_bases":
            # Handle multi-step dissociation for bases
            Kb_values = []
            for key in sorted(compound_data["K"].keys()):
                if key.startswith('Kb'):
                    Kb_values.append(compound_data["K"][key])
            
            OH_minus_total = 0
            for Kb in Kb_values:
                OH_minus = math.sqrt(Kb * concentration)
                OH_minus_total += OH_minus
            
            pOH = -math.log10(OH_minus_total)
            pH = 14 - pOH

        else:
            return jsonify({"error": "Invalid compound category"}), 400

        # Round results to 2 decimal places
        return jsonify({
            "concentration": round(concentration, 4),
            "pH": round(pH, 2),
            "pOH": round(pOH, 2)
        })

    except ValueError as e:
        return jsonify({"error": "Invalid input values"}), 400
    except Exception as e:
        return jsonify({"error": str(e)}), 400

if __name__ == "__main__":
    app.run(debug=True)
