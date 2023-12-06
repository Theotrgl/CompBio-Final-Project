from flask import Flask, request, render_template

app = Flask(__name__)

@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        # Handle the form submission here
        submitted_data = request.form['input_tag']  # Change 'input_name' to the actual name of your input field
        # Process the submitted data or perform any required actions
        print("Submitted data:", submitted_data)  # You can perform any action with the submitted data here
        return "Form submitted successfully!"  

    return render_template('new-1a1o_html_file.html')  

if __name__ == '__main__':
    app.run(debug=True)