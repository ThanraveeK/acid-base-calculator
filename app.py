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
        return jsonify(matches[:10])
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

        # Handle H+ concentration mode
        if calculation_mode == "h-concentration":
            h_conc = float(request.form.get("h_concentration", 0))
            if h_conc <= 0:
                return jsonify({"error": "H+ concentration must be greater than 0"}), 400
            pH = -math.log10(h_conc)
            pOH = 14 - pH
            return jsonify({
                "concentration": h_conc,
                "pH": round(pH, 5),
                "pOH": round(pOH, 5)
            })

        # Handle OH- concentration mode
        elif calculation_mode == "oh-concentration":
            oh_conc = float(request.form.get("oh_concentration", 0))
            if oh_conc <= 0:
                return jsonify({"error": "OH- concentration must be greater than 0"}), 400
            pOH = -math.log10(oh_conc)
            pH = 14 - pOH
            return jsonify({
                "concentration": oh_conc,
                "pH": round(pH, 5),
                "pOH": round(pOH, 5)
            })

        # Handle pH input mode
        elif calculation_mode == "ph-input":
            pH = float(request.form.get("ph_value", 0))
            pOH = 14 - pH
            h_conc = 10 ** (-pH)
            return jsonify({
                "concentration": round(h_conc, 10),
                "pH": round(pH, 5),
                "pOH": round(pOH, 5)
            })

        # Handle pOH input mode
        elif calculation_mode == "poh-input":
            pOH = float(request.form.get("poh_value", 0))
            pH = 14 - pOH
            oh_conc = 10 ** (-pOH)
            return jsonify({
                "concentration": round(oh_conc, 10),
                "pH": round(pH, 5),
                "pOH": round(pOH, 5)
            })

        # Calculate concentration
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
                    return jsonify({"error": "Molecular weight not available"}), 400
                
                moles = mass / molecular_weight
                concentration = moles / volume

        else:
            return jsonify({"error": "Invalid calculation mode"}), 400

        # Calculate pH and pOH
        if compound_category == "strong_acids":
            pH = -math.log10(concentration)
            if concentration == 1:
                pH = 0
            pOH = 14 - pH
            

        elif compound_category == "strong_bases":
            pOH = -math.log10(concentration)
            if concentration == 1:
                pOH = 0
            pH = 14 - pOH

        elif compound_category == "weak_acids":
            Ka = compound_data["K"]
            a = 1
            b = Ka
            c = -Ka * concentration
            H_plus = (-b + math.sqrt(b**2 - 4*a*c)) / (2*a)
            pH = -math.log10(H_plus)
            pOH = 14 - pH

        elif compound_category == "weak_bases":
            Kb = compound_data["K"]
            a = 1
            b = Kb
            c = -Kb * concentration
            OH_minus = (-b + math.sqrt(b**2 - 4*a*c)) / (2*a)
            pOH = -math.log10(OH_minus)
            pH = 14 - pOH

        elif compound_category == "polyprotic_acids":
            K_values = compound_data["K"]
            H_plus_total = 0
            remaining_conc = concentration
            
            for key in sorted(K_values.keys()):
                if key.startswith(('Ka', 'step_')):
                    Ka = K_values[key]
                    H_plus = (-Ka + math.sqrt(Ka*remaining_conc)) 
                    H_plus_total += H_plus
                    remaining_conc -= H_plus
            
            pH = -math.log10(H_plus_total)
            pOH = 14 - pH

        elif compound_category == "polyprotic_bases":
            K_values = compound_data["K"]
            OH_minus_total = 0
            remaining_conc = concentration
            
            for key in sorted(K_values.keys()):
                if key.startswith('Kb'):
                    Kb = K_values[key]
                    OH_minus = (-Kb + math.sqrt(Kb * remaining_conc)) 
                    OH_minus_total += OH_minus
                    remaining_conc -= OH_minus
            
            pOH = -math.log10(OH_minus_total)
            pH = 14 - pOH

        return jsonify({
            "concentration": round(concentration, 5),
            "pH": round(pH, 5),
            "pOH": round(pOH, 5)
        })

    except ValueError as e:
        return jsonify({"error": "Invalid input values"}), 400
    except Exception as e:
        return jsonify({"error": str(e)}), 400

if __name__ == "__main__":
    app.run(debug=True)
