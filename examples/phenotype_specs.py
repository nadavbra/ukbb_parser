import pandas as pd

specs = [
    {
        'name': 'Height',
        'source': 'field',
        'field_id': 50,
        'field_type': 'continuous',
    },
    {
        'name': 'BMI',
        'source': 'field',
        'field_id': 21001,
        'field_type': 'continuous',
    },
    {
        'name': 'Waist circumference',
        'source': 'field',
        'field_id': 48,
        'field_type': 'continuous',
    },
    {
        'name': 'Hip circumference',
        'source': 'field',
        'field_id': 49,
        'field_type': 'continuous',
    },
    {
        'name': 'Diastolic blood pressure',
        'source': 'field',
        'field_id': 4079,
        'field_type': 'continuous',
    },
    {
        'name': 'Systolic blood pressure',
        'source': 'field',
        'field_id': 4080,
        'field_type': 'continuous',
    },
    {
        'name': 'Menarche (age at onset)',
        'source': 'field',
        'field_id': 2714,
        'field_type': 'continuous',
        'sex_filter': 'F',
    },
    {
        'name': 'Menopause (age at onset)',
        'source': 'field',
        'field_id': 3581,
        'field_type': 'continuous',
        'sex_filter': 'F',
    },
    {
        'name': 'Intraocular pressure',
        'source': 'aggregation',
        'subspecs': [
            {
                'name': 'Intraocular pressure (right)',
                'source': 'field',
                'field_id': 5254,
                'field_type': 'continuous',
            },
            {
                'name': 'Intraocular pressure (left)',
                'source': 'field',
                'field_id': 5262,
                'field_type': 'continuous',
            },
        ],
        'aggregation_function': lambda right, left: pd.concat([right, left], axis = 1).max(axis = 1), # We take the maximum between the two fields.
    },
    {
        'name': 'Hand grip strength',
        'source': 'aggregation',
        'subspecs': [
            {
                'name': 'Hand grip strength (left)',
                'source': 'field',
                'field_id': 46,
                'field_type': 'continuous',
            },
            {
                'name': 'Hand grip strength (right)',
                'source': 'field',
                'field_id': 47,
                'field_type': 'continuous',
            },
        ],
        'aggregation_function': lambda left, right: pd.concat([left, right], axis = 1).max(axis = 1), # We take the maximum between the two fields.
    },
    {
        'name': 'Male-pattern baldness',
        'source': 'field',
        'field_id': 2395,
        'field_type': 'set',
        'one_hot_encoding': True,
        'sex_filter': 'M',
    },
    {
        'name': 'Platelet count',
        'source': 'field',
        'field_id': 30080,
        'field_type': 'continuous',
    },
    {
        'name': 'Monocyte count',
        'source': 'field',
        'field_id': 30130,
        'field_type': 'continuous',
    },
    {
        'name': 'Red blood cell count',
        'source': 'field',
        'field_id': 30010,
        'field_type': 'continuous',
    },
    {
        'name': 'White blood cell count',
        'source': 'field',
        'field_id': 30000,
        'field_type': 'continuous',
    },
    {
        'name': 'High light scatter reticulocyte count',
        'source': 'field',
        'field_id': 30300,
        'field_type': 'continuous',
    },
    {
        'name': 'Eosinophil counts',
        'source': 'field',
        'field_id': 30150,
        'field_type': 'continuous',
    },
    {
        'name': 'Reticulocyte count',
        'source': 'field',
        'field_id': 30250,
        'field_type': 'continuous',
    },
    {
        'name': 'Lymphocyte counts',
        'source': 'field',
        'field_id': 30120,
        'field_type': 'continuous',
    },
    {
        'name': 'Mean platelet volume',
        'source': 'field',
        'field_id': 30100,
        'field_type': 'continuous',
    },
    {
        'name': 'Mean corpuscular volume',
        'source': 'field',
        'field_id': 30040,
        'field_type': 'continuous',
    },
    {
        'name': 'Mean corpuscular hemoglobin',
        'source': 'field',
        'field_id': 30050,
        'field_type': 'continuous',
    },
    {
        'name': 'Neutrophil count',
        'source': 'field',
        'field_id': 30140,
        'field_type': 'continuous',
    },
    {
        'name': 'Red cell distribution width',
        'source': 'field',
        'field_id': 30070,
        'field_type': 'continuous',
    },
    {
        'name': 'Platelet distribution width',
        'source': 'field',
        'field_id': 30110,
        'field_type': 'continuous',
    },
    {
        'name': 'High light scatter reticulocyte percentage of red cells',
        'source': 'field',
        'field_id': 30290,
        'field_type': 'continuous',
    },
    {
        'name': 'Breast cancer',
        'source': 'ICD-10',
        'codings': ['C50'],
        'sex_filter': 'F',
    },
    {
        'name': 'Epithelial ovarian cancer',
        'source': 'ICD-10',
        'codings': ['C56'],
        'sex_filter': 'F',
    },
    {
        'name': 'Prostate cancer',
        'source': 'ICD-10',
        'codings': ['C61'],
        'sex_filter': 'M',
    },
    {
        'name': 'Colorectal cancer',
        'source': 'ICD-10',
        'codings': ['C18'],
    },
    {
        'name': 'Lung cancer',
        'source': 'ICD-10',
        'codings': ['C34'],
    },
    {
        'name': 'Chronic lymphocytic leukemia',
        'source': 'ICD-10',
        'codings': ['C91'],
    },
    {
        'name': 'Pancreatic cancer',
        'source': 'ICD-10',
        'codings': ['C25'],
    },
    {
        'name': 'Melanoma',
        'source': 'ICD-10',
        'codings': ['C43'],
    },
    {
        'name': 'Schizophrenia',
        'source': 'ICD-10',
        'codings': ['F20'],
    },
    {
        'name': 'Bipolar disorder',
        'source': 'ICD-10',
        'codings': ['F31'],
    },
    {
        'name': 'Major depressive disorder',
        'source': 'ICD-10',
        'codings': ['F33'],
    },
    {
        'name': 'Parkinson\'s disease',
        'source': 'ICD-10',
        'codings': ['G20'],
    },
    {
        'name': 'Stroke',
        'source': 'ICD-10',
        'codings': ['I63'],
    },
    {
        'name': 'Hypertension',
        'source': 'ICD-10',
        'codings': ['I10'],
    },
    {
        'name': 'Sudden cardiac arrest',
        'source': 'ICD-10',
        'codings': ['I46'],
    },
    {
        'name': 'Type 1 diabetes',
        'source': 'ICD-10',
        'codings': ['E10'],
    },
    {
        'name': 'Type 2 diabetes',
        'source': 'ICD-10',
        'codings': ['E11'],
    },
    {
        'name': 'Systemic sclerosis',
        'source': 'ICD-10',
        'codings': ['M34'],
    },
    {
        'name': 'Multiple sclerosis',
        'source': 'ICD-10',
        'codings': ['G35'],
    },
    {
        'name': 'Systemic lupus erythematosus',
        'source': 'ICD-10',
        'codings': ['M32'],
    },
    {
        'name': 'Rheumatoid arthritis',
        'source': 'ICD-10',
        'codings': ['M05', 'M06'],
    },
    {
        'name': 'Asthma',
        'source': 'ICD-10',
        'codings': ['J45'],
    },
    {
        'name': 'Crohn\'s and colitis',
        'source': 'ICD-10',
        'codings': ['K50', 'K51'],
    },
]