# Additional tests to add into test_obj_defaults.py.  Not yet
# ready for prime-time.



    def test_Qtcm_object_field_default_attributes_on_init_no_extras(self):
        """Test Qtcm object no extra field attributes on init.
        """
        model = Qtcm()
        self.failUnless( not hasattr(model, 'viscU') )
        self.failUnless( not hasattr(model, 'viscT') )
        self.failUnless( not hasattr(model, 'viscQ') )


    def test_Qtcm_object_compiled_form_attr(self):
        """Test Qtcm object compiled_form attribute and keyword init.
        """
        inputs = {}
        inputs['dt'] = 600.
        inputs['eps_c'] = 0.15
        model = Qtcm(**inputs)
        self.failUnlessEqual( model.compiled_form, 'full' )
        del model

        inputs['compiled_form'] = 'full'
        model = Qtcm(**inputs)
        self.failUnlessEqual( model.compiled_form, 'full' )
        del model

        inputs = {}
        inputs['compiled_form'] = 'oops'
        self.failUnlessRaises(ValueError, Qtcm, **inputs)

        inputs = {}
        model = Qtcm(**inputs)
        model.compiled_form = 'hello'
        self.failUnlessRaises(ValueError, model.run)
        del model


    def test_Qtcm_object_get_items_from_compiled(self):
        """Test Qtcm get_qtcm_item method.
        """
        inputs = {}
        inputs['dt'] = 120.
        inputs['eps_c'] = 0.15
        model = Qtcm(**inputs)
        self.failUnless( self.N.allclose(model.get_qtcm_item('dt'), 120.) )
        self.failUnless( self.N.allclose(model.get_qtcm_item('eps_c'), 0.15) )
        del model

        model = Qtcm()
        self.failUnless( self.N.allclose(model.get_qtcm_item('dt'), 1200.) )
        self.failUnless( \
           self.N.allclose(model.get_qtcm_item('eps_c'), 0.13888889E-03) )
        self.failUnlessEqual( model.get_qtcm_item('title') \
                            , 'QTCM default title' )
        self.failUnlessEqual( model.get_qtcm_item('bnddir'), '../bnddata' )
        self.failUnlessEqual( model.get_qtcm_item('SSTdir') \
                            , '../bnddata/SST_Reynolds' )
        self.failUnlessEqual( model.get_qtcm_item('outdir') \
                            , '../proc/qtcm_output' )
        self.failUnlessEqual( model.get_qtcm_item('runname'), 'runname' )
        self.failUnlessEqual( model.get_qtcm_item('landon'), 1 )
        self.failUnlessEqual( model.get_qtcm_item('SSTmode'), 'seasonal' )
        self.failUnlessEqual( model.get_qtcm_item('year0'), 0 )
        self.failUnlessEqual( model.get_qtcm_item('month0'), -1 )
        self.failUnlessEqual( model.get_qtcm_item('day0'), -1 )
        self.failUnlessEqual( model.get_qtcm_item('lastday'), 365 )
        self.failUnlessEqual( model.get_qtcm_item('interval'), 1 )
        self.failUnlessEqual( model.get_qtcm_item('noout'), 0 )
        self.failUnlessEqual( model.get_qtcm_item('nooutr'), 0 )
        self.failUnlessEqual( model.get_qtcm_item('ntout'), -30 )
        self.failUnlessEqual( model.get_qtcm_item('ntouti'), 0 )
        self.failUnlessEqual( model.get_qtcm_item('ntoutr'), 0 )
        self.failUnlessEqual( model.get_qtcm_item('mrestart'), 1 )
        self.failUnlessEqual( model.get_qtcm_item('mt0'), 1 )
        self.failUnless( self.N.allclose(model.get_qtcm_item('ziml'), 500.) )
        self.failUnless( self.N.allclose(model.get_qtcm_item('weml'), 0.01) )
        self.failUnless( self.N.allclose(model.get_qtcm_item('VVsmin'), 4.5) )
        self.failUnlessEqual( model.get_qtcm_item('arr1name'), '?' )
        self.failUnlessEqual( model.get_qtcm_item('arr2name'), '?' )
        self.failUnlessEqual( model.get_qtcm_item('arr3name'), '?' )
        self.failUnlessEqual( model.get_qtcm_item('arr4name'), '?' )
        self.failUnlessEqual( model.get_qtcm_item('arr5name'), '?' )
        self.failUnlessEqual( model.get_qtcm_item('arr6name'), '?' )
        self.failUnlessEqual( model.get_qtcm_item('arr7name'), '?' )
        self.failUnlessEqual( model.get_qtcm_item('arr8name'), '?' )
        self.failUnless( self.N.allclose(model.get_qtcm_item('viscxu0'), 7e5) )
        self.failUnless( self.N.allclose(model.get_qtcm_item('viscyu0'), 7e5) )
        self.failUnless( self.N.allclose(model.get_qtcm_item('visc4x'), 7e5) )
        self.failUnless( self.N.allclose(model.get_qtcm_item('visc4y'), 7e5) )
        self.failUnless( self.N.allclose(model.get_qtcm_item('viscxu1'), 7e5) )
        self.failUnless( self.N.allclose(model.get_qtcm_item('viscyu1'), 7e5) )
        self.failUnless( self.N.allclose(model.get_qtcm_item('viscxT'), 12e5) )
        self.failUnless( self.N.allclose(model.get_qtcm_item('viscyT'), 12e5) )
        self.failUnless( self.N.allclose(model.get_qtcm_item('viscxq'), 12e5) )
        self.failUnless( self.N.allclose(model.get_qtcm_item('viscyq'), 12e5) )
        del model


    def test_Qtcm_object_set_and_get_items_from_compiled(self):
        """Test Qtcm set_qtcm_item and get_qtcm_item methods together.
        """
        model = Qtcm()
        model.set_qtcm_item('dt', 2400.)
        model.set_qtcm_item('bnddir', 'ooga booga')
        model.set_qtcm_item('ntout', 120)

        self.failUnless( self.N.allclose(model.get_qtcm_item('dt'), 2400.) )
        self.failUnlessEqual( model.get_qtcm_item('bnddir'), 'ooga booga' )
        self.failUnlessEqual( model.get_qtcm_item('ntout'), 120 )

        self.failUnlessRaises(TypeError, model.set_qtcm_item, 'dt', 400)
        self.failUnlessRaises(TypeError, model.set_qtcm_item, 'dt', 'hi')
        self.failUnlessRaises(TypeError, model.set_qtcm_item, 'bnddir', 20)
        self.failUnlessRaises(TypeError, model.set_qtcm_item, 'bnddir', 20.)
        self.failUnlessRaises(TypeError, model.set_qtcm_item, 'ntout', 30.)
        self.failUnlessRaises(TypeError, model.set_qtcm_item, 'ntout', 'hi')


    def test_Qtcm_object_sync_all_py_values_to_qtcm_items(self):
        """Test Qtcm set_qtcm_item and get_qtcm_item methods together.
        """
        #- Check that using the private method _set_qtcm_item_in_model
        #  does not change the value of the attribute from the default:

        model = Qtcm()
        model._set_qtcm_item_in_model('dt', 240.420)
        model._set_qtcm_item_in_model('bnddir', 'ooga booga')
        model._set_qtcm_item_in_model('ntout', 12320)
        self.failUnless( self.N.allclose(model.dt.value, Field('dt').value) )
        self.failUnlessEqual( model.bnddir.value, Field('bnddir').value )
        self.failUnlessEqual( model.ntout.value, Field('ntout').value )


    def test_Qtcm_object_only_one_instance(self):
        """Test Qtcm only can have one instance.
        """
        model = Qtcm()
        self.failUnlessRaises(RuntimeError, Qtcm)
        del model

        model = Qtcm(compiled_form='full')
        self.failUnlessRaises(RuntimeError, Qtcm)
        del model


    def test_so_not_changed_permanently_between_objects(self):
        """Test .so objects not changed permanently between instantiations.
        """
        model = Qtcm()
        model._set_qtcm_item_in_model('dt', 240.420)
        model._set_qtcm_item_in_model('bnddir', 'ooga booga')
        model._set_qtcm_item_in_model('ntout', 12320)
        m1_dt = copy.copy(model.get_qtcm_item('dt'))
        m1_bnddir = copy.copy(model.get_qtcm_item('bnddir'))
        m1_ntout = copy.copy(model.get_qtcm_item('ntout'))
        del model

        model = Qtcm()
        self.failIfAlmostEqual( m1_dt, model.get_qtcm_item('dt') )
        self.failIfEqual( m1_bnddir, model.get_qtcm_item('bnddir') )
        self.failIfEqual( m1_ntout, model.get_qtcm_item('ntout') )
