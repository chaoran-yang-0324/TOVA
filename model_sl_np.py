#!/usr/bin/env python
# coding: utf-8

# In[2]:


# pip install streamlit


# In[1]:


import streamlit as st
import numpy as np

def long_function(input1, input2, input3):
    # Example of your function using numpy
    return np.array([input1 * 2, input2 + 3, input3 - 1])

st.title("Interactive Function App")
st.write("Enter the values for input1, input2, and input3:")

input1 = st.number_input("Input 1", value=0.0)
input2 = st.number_input("Input 2", value=0.0)
input3 = st.number_input("Input 3", value=0.0)

if st.button("Calculate"):
    output = long_function(input1, input2, input3)
    st.write(f"Output: {output}")

