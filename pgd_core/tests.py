"""
This file demonstrates writing tests using the unittest module. These will pass
when you run "manage.py test".

Replace this with more appropriate tests for your application.
"""

from django.test import TestCase , Client
from django.core.urlresolvers import reverse
from django.core import mail


class RegistrationTestCase(TestCase):


	fixtures    	  = ['users']
	default_url 	  = "http://testserver"
	test_email		  = {'email' : 'email@example.org'}
	test_ceredentials = {'username':'test1', 'password':'hey'}
	
	register_details  = {
	'username' : 'vamos', 
	'email' : 'email@example.org', 
	'password1' : 'hola',
	'password2' : 'hola'}

	test_new_password = {'old_password': 'hey','new_password1': 'heythere',
						'new_password2' : 'heythere'}

	def test_register(self) :
		test_client = Client()
		response = test_client.post(reverse('registration_register'), self.register_details, follow=True)
		self.assertEqual(response.status_code, 200)
		self.assertEqual(response.redirect_chain[-1][0], self.default_url+'/accounts/register/complete/')
		self.assertEqual(len(mail.outbox), 1)
		self.assertEqual(mail.outbox[0].to[0], self.test_email['email'])	

	def test_reset_password(self):
		response = self.client.post(reverse('auth_password_reset'), self.test_email)
		self.assertEqual(response.status_code, 302)
		self.assertEqual(response['Location'], self.default_url+reverse('auth_password_reset_done'))

	def test_change_password(self):
		test_client = Client()
		get_response = test_client.get(reverse('auth_password_change'))
		self.assertEqual(get_response.status_code, 302)
		post_response_1 = test_client.post(get_response['Location'], self.test_ceredentials, follow=True)
		self.assertEqual(post_response_1.status_code, 200)
		post_response_2 = test_client.post(post_response_1.redirect_chain[-1][0], self.test_new_password, follow=True)
		self.assertEqual(post_response_2.status_code, 200)
		self.assertEqual(post_response_2.redirect_chain[-1][0], self.default_url+reverse('auth_password_change_done'))